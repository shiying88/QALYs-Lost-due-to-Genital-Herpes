import re
import pandas as pd


# get x from string represents interval from x to (x+1)
def get_initial_age(cell):
    split_list = re.split('[– ]', cell)
    lower_age = float(split_list[0])
    # upper_age = float(split_list[-1])
    return lower_age


# convert string (with ',' as delimiter) to int
def convert_to_int(num_string):
    num = int(num_string.replace(',', ''))
    return num


# get number surviving at age x, convert number with comma format to int format
def get_number_surviving(df, age):
    row = df.loc[df['Start age'] == age]
    num_string = list(row['Number surviving to age x'])[0]
    num_int = convert_to_int(num_string=num_string)
    return num_int


# read dataset, preprocess data
def get_clean_life_table(path=None, df_filename=None):
    if path is None:
        df = df_filename
    else:
        df = pd.read_csv(path)
    # drop rows with missing values
    df = df.dropna()
    # replace header with the top row
    new_header = df.iloc[0]  # grab the first row for the header
    df = df[1:]  # take the data less the header row
    df.columns = new_header  # set the header row as the df header
    # add new column for starting age and length of interval
    df['Start age'] = df['Age (years)'].apply(get_initial_age)
    interval_length = [1] * len(df)
    df['Interval length'] = interval_length
    # column names of the outcome dataframe:
    # 'Age (years)', 'Probability of dying between ages x and x + 1',
    # 'Number surviving to age x', 'Number dying between ages x and x + 1',
    # 'Person–years lived between ages x and x + 1',
    # 'Total number of person–years lived above age x',
    # 'Expectation of life at age x',
    # 'Start age', 'Interval length'
    return df
