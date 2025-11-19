import requests
import numpy as np
import time
from pathlib import Path
from io import StringIO
import pandas as pd

path = Path(__file__).parent.absolute()

def server_request(user_params,url_mod):
    """This function requests data from the API.

    Args:
        user_params (dictionary): The parameters for the request, will usually be filter and sort.
        url_mod (string): Where to request from, e.g. launches, objects, vehicles etc. 

    Returns:
        response: The response from the server (deal with errors in main function).
    """
    
    #Define the URL for the data, and the authorisation token (get a new token here if needed https://discosweb.esoc.esa.int/tokens).
    api_url = 'https://discosweb.esoc.esa.int/api'
    token = np.loadtxt(str(path)+"/DISCOSweb_token.txt", dtype=str)

    response = requests.get(
            #Set up the url name.
            f'{api_url}{url_mod}',
            #Set up authorization.
            headers={
                'Authorization': f'Bearer {token}',
                'DiscosWeb-Api-Version': '2',
            },
            #Use the user-specified parameters.
            params=user_params
        ) 
    return response

def wait_function(message,wait_time):
    """This function outputs a message to the screen telling the user the wait time.

    Args:
        message (string): An additional message to print to the screen about the current process.
        wait_time (int): How many seconds are left befor emore requests can be made.
    """
    print("\n")        
    for i in range(wait_time,0,-1):
        #The 'end="\r", flush=True' means that the line will replace itself instead of printing 60 lines.
        #This is useful to clean up the terminal, but may look weird sometimes before/after waiting.
        print(f"{message}Too many requests, retrying in {str(i)}s.", end="\r", flush=True)
        time.sleep(1)
    print("\n")
    
def response_error_handler(response,message):
    if response.status_code == 400:
        print("Client Error")
    elif response.status_code == 429:
        wait_time=(int(response.headers["X-Ratelimit-Reset"])-int(time.time()))+1
        wait_function(message,wait_time)#

def scrape_jsr(url,session):
    # Function to web scrape data and convert to a pandas DataFrame.
    response = session.get(url)
    if response.status_code == 200:
        # Convert the content to a file-like object for pandas
        tsv_data = StringIO(response.text)
        # Load the data into a pandas DataFrame
        df = pd.read_csv(tsv_data, delimiter="\t", dtype=object, low_memory=False)
    else:
        raise ImportError(f"Failed to fetch from JSR URL: {url}", response.status_code)
    return df