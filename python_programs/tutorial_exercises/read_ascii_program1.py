fileName = '/Users/md875/data/tutorial_data/smi_order.txt'
fileData = open(fileName, 'r')

all_lines = fileData.readlines()  #reads entire ascii file as one very long stream of continuous ascii chars as a list
fileData.close()
print all_lines

all_lines = ''.join(all_lines)  #converts the list to string with elements separated by ''
print all_lines

fileName1 = '/Users/md875/data/tutorial_data/smi_order_copy.txt'
fileData1 = open(fileName1, 'w')

fileData1.write('The quick brown fox')
fileData1.write(' jumped over the \n')
fileData1.write('lazy dog.')
fileData1.close()


