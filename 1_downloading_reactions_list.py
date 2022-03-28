## Libraries
import requests


## Getting reactions list
r = requests.post("http://rest.kegg.jp/list/reaction")

## Saving reactions list
with open('data/kegg_reactions.txt', mode='w') as myFile:
    myFile.write(r.text)
