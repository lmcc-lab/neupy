from scipy.linalg import sqrtm
import numpy as np
import parameters
import time

params = parameters.read_params()

# Gather parameters from paratmeters file
ALPHA = params.ALPHA
BETA = params.BETA
KAPPA = params.KAPPA



class UT():

    def __init__(self, x, eta, P, Q):
        """
        Unscented tranform class for propegating uncertainties through non-linear equations

        :param x: initial state vector [r, v] (n x 1)
        :param eta: initial state noise vector (q x 1)
        :param P: initial covariance matrix (n x n)
        :param Q: initial processs noise covariance matrix (q x q)

        :returns    self.xa: Augmented state vector [x.T, eta.T].T
                    self.Pa: Augmented covariance matrix [[P, 0],
                                                          [0,Q]]
                    self.na: Number of dimensions in the augmented state vector (n+q)
                    self.n: Number of dimensions in the state vector n
                    self.q: Number of dimensions in the noise vector q
        """
        self.x = x
        self.eta = eta
        self.P = P
        self.Q = Q

        self.xa = np.vstack((self.x, self.eta))
        PaUp = np.hstack((self.P, np.zeros((self.P.shape[0], self.Q.shape[1]))))
        PaDw = np.hstack((np.zeros((self.Q.shape[0], self.P.shape[1])), self.Q))
        self.Pa = np.vstack((PaUp, PaDw))
        self.na = len(self.xa[:,0])
        self.n = len(self.x[:,0])
        self.q = len(self.eta[:,0])

    def weights(self):
        """
        Calculate weights and gamma.

        :return: wm mean weights
                 wc covariance weights
                 gamma scaling factor.
        """


        lam = ALPHA**2*(self.na+KAPPA)-self.na

        wm = [lam/(self.na+lam)]

        wc = [lam/(self.na+lam+1-ALPHA**2+BETA)]

        for i in range(2*self.na):
            wm.append(1.0/(2*(self.na+lam)))
            wc.append(1.0/(2*(self.na+lam)))

        gamma = np.sqrt(self.na+lam)

        wm = np.array(wm)
        wc = np.array(wc)

        return wm, wc, gamma

    def calc_sigma_points(self, gamma):
        """
        Calculate the sigma points used for the propegation of the uncertainties.
        This generates the sigma matrix, which is seperated into a state sigma
        matrix and a process noise sigma matrix.

        :param gamma: sqrt(n+lam) from weights function

        :returns    self.x_sigma: state sigma matrix
                    self.w_sigma: noise sigma matrix
                    self.a_sigma: augmented sigma matrix
        """


        self.a_sigma = self.xa                                                                  # Initialise self.a_sigma (augmented sigma points)

        L = sqrtm(self.Pa)                                                                      # Calculate the square root of the augmented covariance matrix

        for i in range(self.na):
            self.a_sigma = np.hstack((self.a_sigma, self.xa + gamma * L[:, i:i+1]))             # horizontally stack the first na augmented sigma points.

        for i in range(self.na):
            self.a_sigma = np.hstack((self.a_sigma, self.xa - gamma * L[:, i:i+1]))             # horizontally stack the next na+1 to 2na augmented sigma points

        self.x_sigma = self.a_sigma[:self.n, :]                                                 # Seperate the augmented sigma matrix into its state sigma matrix component.
        self.w_sigma = self.a_sigma[self.n:, :]                                                 # Seperate the augmented sigma matrix into its noise sigma matrix component.

    def propegate_sigma_points(self, motion_model):
        """
        Propegate sigma points through non-linear motion model

        :param motion_model: function where the motion model is held
        """

        for i in range(self.x_sigma.shape[1]):
            self.x_sigma[:,i:i+1] = motion_model(self.x_sigma[:,i:i+1], self.w_sigma[:,i:i+1])  # x_sigma are the states of the sigmas, while w_sigma are the uncertainties.

    def prediction_step(self, wm, wc, motion_model):
        """
        The prediction step propegates the sigma points at time step k to time step k+1,
        calculates the new mean and propegates the covariance matrix.

        :param wm: (1,) array - weights of means
        :param wc: (1,) array - weights of covariance matrices
        :param motion_model: function - motion model for propegation

        """


        self.propegate_sigma_points(motion_model)                                               # Propegate sigma points through the motion model from time step k to k+1

        self.x = (wm @ self.x_sigma.T).T                                                        # Calculate the new mean state vector

        self.P = np.zeros((self.n, self.n))                                                     # Initialise the covariance matrix

        d = self.x_sigma.T - self.x                                                             # Calculate the distance between the sigma points and the new mean

        d = d.T                                                                                 # Transpose for correct orientation
        for i in range(2*self.na):
            self.P = self.P + wc[i]* d[:, i:i+1] @ d[:, i:i+1].T                                # Calculate the covariance matrix performing the weighted sum of all of the distance matrices.

        self.x = np.array([self.x]).T                                                           # Transpose the state matrix back to its correct orientation.










