import csv 
import glob

files = glob.glob("../data/samples/*")


for file in files:  
    name = file.split('-60k.csv')[0].split('/')[-1]
    print(name)
    with open(file, 'rb') as f: 
        vec = f.readlines()
    aux = []
    for i in range(2, len(vec)):
        aux.append(int(vec[i]) - int(vec[i - 1]))
    
    print(len(aux))
    for i in range(2, len(aux) -1):
        writer = open(f"scatter/{name}.txt", 'a')
        writer.write(f"{aux[i]} {aux[i - 1]}\n")
        writer.close()
