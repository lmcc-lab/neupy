import pandas as pd

test = {"col1": [1, 2, 3], "col2": [3, 4, 4]}

df = pd.DataFrame(test)

for row in df[::-1]:
    print(row)