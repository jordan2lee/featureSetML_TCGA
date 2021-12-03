import json

#####
f_out = 'src/modelID_performance2importance_ALLCOHORTS.json'
#####
data = {}
with open('src/conversions/subSCOPE.json', 'r') as subscope:
    data.update(json.load(subscope))
with open('src/conversions/skgrid.json', 'r') as subscope:
    data.update(json.load(subscope))
with open('src/conversions/CF.json', 'r') as subscope:
    data.update(json.load(subscope))
with open('src/conversions/aklimate.json', 'r') as subscope:
    data.update(json.load(subscope))
with open('src/conversions/jadbio.json', 'r') as subscope:
    data.update(json.load(subscope))

# Output conversion keys
with open(f_out, 'w') as out:
    out.write(json.dumps(data))
    out.write('\n')
