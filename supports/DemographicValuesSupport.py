import re
import pandas as pd


####################
# Helper functions #
####################
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


###########################
# disutility due to aging #
###########################
# data source: EQ-5D background disutility estimation
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2634296/
def get_bg_utility(age):
    """ get background utility value based on age """
    bg_u = None  # initial background utility
    if age < 18:
        bg_u = 1
    if 18 <= age < 30:
        bg_u = 0.922
    elif 30 <= age < 40:
        bg_u = 0.901
    elif 40 <= age < 50:
        bg_u = 0.871
    elif 50 <= age < 60:
        bg_u = 0.842
    elif 60 <= age < 70:
        bg_u = 0.823
    elif 70 <= age < 80:
        bg_u = 0.790
    elif age >= 80:
        bg_u = 0.736
    return bg_u


def get_adjusted_disu(current_age, disu, out='disu'):
    """ get adjusted disutility based on disease disutility and background disutility due to aging
    :param current_age: current age of people
    :param disu: represents disutility of disease
    :param out: 'disu': we want to get adjusted disutility; 'utl': we want to get adjusted utility
    """
    bg_u = get_bg_utility(age=current_age)
    adjusted_disu = bg_u * disu
    adjusted_utl = bg_u - adjusted_disu
    if out == 'utl':
        return adjusted_utl
    else:
        return adjusted_disu


##########################
#       life table       #
##########################
# read life table for male and female, 2017
# link: https://www.cdc.gov/nchs/data/nvsr/nvsr67/nvsr67_07-508.pdf
df_life_table_male = get_clean_life_table(
    path='source/data/US_man_lifetable2017_CDC.csv')
df_life_table_female = get_clean_life_table(
    path='source/data/US_women_lifetable2017_CDC.csv')
df_life_table = get_clean_life_table(
    path='source/data/US_lifetable2015_CDC.csv')


def get_conditional_survival_rate(younger_age, older_age, sex=None):
    younger_age = round(younger_age)
    older_age = round(older_age)
    if sex == 'female':
        df = df_life_table_female
    elif sex == 'male':
        df = df_life_table_male
    elif sex is None:
        df = df_life_table
    else:
        raise ValueError('wrong sex type')
    # calculate conditional surviving rate at age_old giving people being alive at age_young
    num_alive_young = get_number_surviving(df=df, age=younger_age)
    num_alive_old = get_number_surviving(df=df, age=older_age)
    conditional_prob = num_alive_old / num_alive_young
    return round(conditional_prob, 4)


# the target group aged between 18-49
def create_new_life_table_for_QALE(df, seq_disu, discount=0.03, excess_mort=None):
    """ create a new DataFrame for quality adjusted life expectancy calculation
    :param df: DataFrame
    :param seq_disu: disutility of sequelae
    :param discount: rate for discounting
    :param excess_mort: excess mortality
    :return (DataFrame) with columns of:
    age, death_rate, alive, bg_u, adj_utl, uw_nd, uw_d
    """
    df_new = pd.DataFrame()
    df_new['age'] = df['Start age']
    df_new['death_rate'] = df['Probability of dying between ages x and x + 1']
    if excess_mort is not None:
        df_new['death_rate'] = df_new.apply(lambda row: (row['death_rate'] + excess_mort), axis=1)
    # add proportion of people alive
    current_alive = 1
    new_col = []
    for ind in df_new.index:
        new_col.append(current_alive)
        current_alive = current_alive * (1 - float(df_new['death_rate'][ind]))
    df_new['alive'] = new_col
    # background utility weights based on age
    df_new['bg_utl'] = df_new.apply(lambda row: get_bg_utility(age=row['age']), axis=1)
    # adjusted utility (= bg_utl - bg adjusted seq_disu)
    df_new['adj_utl'] = df_new.apply(lambda row: row['bg_utl'] * (1 - seq_disu), axis=1)
    # utility-weights without discounting
    df_new['uw_nd'] = df_new.apply(lambda row: row['alive'] * row['adj_utl'], axis=1)
    # utility-weights with discounting
    df_new['uw_d'] = df_new.apply(lambda row: (row['uw_nd'] / (1 + discount) ** row['age']), axis=1)

    return df_new


def get_QALE(seq_disu, discount=0.03, excess_mort=None):
    """
    life table approach to estimate life expectancy for neonates
    :param seq_disu: (float) disutility of nor/mil/mod/sev
    :param discount: (float) default None (no discounting)
    :param excess_mort: (float) excess mortality caused by neurological sequelae
    :return: (float) life expectancy
    """
    # read original life table
    df_female = df_life_table_female
    df_male = df_life_table_male
    # sex ratio at birth
    male_ratio = 0.511480215
    female_ratio = 0.488519785
    # get new DataFrame
    df_female_new = create_new_life_table_for_QALE(df=df_female, seq_disu=seq_disu,
                                                   discount=discount, excess_mort=excess_mort)
    df_male_new = create_new_life_table_for_QALE(df=df_male, seq_disu=seq_disu,
                                                 discount=discount, excess_mort=excess_mort)
    print(round(sum(df_male_new['uw_nd']), 1))
    print(round(sum(df_male_new['uw_d']), 1))
    # calculate life expectancy
    d_life_expectancy = sum(df_male_new['uw_d']) * male_ratio + sum(df_female_new['uw_d']) * female_ratio
    nd_life_expectancy = sum(df_male_new['uw_nd']) * male_ratio + sum(df_female_new['uw_nd']) * female_ratio
    return nd_life_expectancy, d_life_expectancy
