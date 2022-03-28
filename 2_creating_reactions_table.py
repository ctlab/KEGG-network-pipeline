## Libraries
import csv
import re
import requests
import time


## Getting reactions list
with open('data/kegg_reactions.txt') as myFile:
    contents = myFile.read()

kegg_reactions_list = re.findall(r'rn:(R\d+)', str(contents))
print("Total amount of KEGG reactions:", len(kegg_reactions_list))


## Getting reactions in EQUATION format
start_time = time.time()

kegg_reactions_data = [["Reaction ID", "Reaction", "Enzyme"]]
count = 0


for kegg_reaction_id in kegg_reactions_list:
    annotation_element = []

    url = "http://rest.kegg.jp/get/" + kegg_reaction_id
    try:
        r = requests.post(url)
    except ConnectionResetError:
        time.sleep(3)
        r = requests.post(url)
    except URLError:
        time.sleep(6)
        r = requests.post(url)
    except:
        time.sleep(10)
        r = requests.post(url)

    equation = re.findall(r'EQUATION[ ]+(.+)', str(r.text))
    enzyme = re.findall(r'ENZYME[ ]+(.+)', str(r.text))

    if enzyme != []:
        kegg_reactions_data.append([kegg_reaction_id, equation[0], enzyme[0]])
    else:
        kegg_reactions_data.append([kegg_reaction_id, equation[0], "-"])

    count += 1

    # TEST USE ONLY
    if count == 10:
        break


# Saving reactions as csv
myFile = open('data/kegg_reactions.csv', 'w', newline='')
with myFile:
    writer = csv.writer(myFile, delimiter=';')
    writer.writerows(kegg_reactions_data)
