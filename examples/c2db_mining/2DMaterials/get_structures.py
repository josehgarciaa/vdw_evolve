import sys
sys.path.append('../')
import lib.parser as c2db_parser


uids = []
with open("2Dmagnets_uids.dat", 'r') as f:   
    while True :
        uid = f.readline().strip()
        uids.append(uid)
        if len(uid)==0:
            break

# Print
print("The number of elements is:", len(uids))
c2db_parser.extract_structure(uids)
