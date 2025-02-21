import requests

# ADS_TOKEN = "Ft57LMFzzfeBH6ByOlXmufZhadIDJ4ay50n4o3mJ"  # Use your actual API key
# BIBCODE =  "2024A%26A...688A..69F"  # Replace with your real bibcode
# headers = {"Authorization": f"Bearer {ADS_TOKEN}"}
#
# url = f"https://api.adsabs.harvard.edu/v1/search/query?q=bibcode:{BIBCODE}&fl=citation_count"
# headers = {"Authorization": f"Bearer {ADS_TOKEN}"}
#
# response = requests.get(url, headers=headers)
# print(response.json())  # Print the full response for debugging
#
#
# # import requests

ADS_TOKEN = "Ft57LMFzzfeBH6ByOlXmufZhadIDJ4ay50n4o3mJ"
LIBRARY_BIBCODE = "2024A%26A...688A..69F"  # Replace with your library's ADS bibcode

def get_citation_count():
    headers = {"Authorization": 'Bearer Ft57LMFzzfeBH6ByOlXmufZhadIDJ4ay50n4o3mJ'}
    url = f"https://api.adsabs.harvard.edu/v1/search/query?q=bibcode:{LIBRARY_BIBCODE}&fl=citation_count"
    # url = f"https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..69F/abstract"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        return data["response"]["docs"][0]["citation_count"]
        # return data["response"]
    return "N/A"

print(get_citation_count())

# with open("docs/citation_count.txt", "w") as f:
#     f.write(str(citation_count))
