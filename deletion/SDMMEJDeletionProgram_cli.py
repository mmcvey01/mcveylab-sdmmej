from EJR_junction_auxfns import *
from EJR_junction_pool import *
from EJR_junction_class import *
import argparse
import pandas as pd
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hi", dest="h_in", help="hifiber results, processed by previous step", required=True)
    parser.add_argument("-del", dest="d_in", help="deletion file", required=True)
    parser.add_argument("-out", dest="dir_path", help="directory_for_outputs", required=True)
    parser.add_argument("-n", dest="nick", help="location of the nick (the 0-based reference base after which nick occurs", required=True)
    parser.add_argument("-r", dest="search_radius", help="number of bases to the left and right of the repair junction to search for SD-MMEJ consistent repeats", required=False, default=30)
    parser.add_argument("-m", dest="min_consist_rep_size", help="minimum size of an SD-MMEJ consistent repeat", required=False, default=4)

    args = parser.parse_args()

    return args.h_in, args.d_in, args.dir_path, int(args.nick), int(args.search_radius), int(args.min_consist_rep_size)

def parse_ref(h_in, nick):
    hifi_input = pd.read_csv(h_in)

    reference = hifi_input.loc[hifi_input['CLASS'] == 'exact']

    if (len(reference) != 1):
        print(len(reference))
        print("ERROR wrong number of reference sequences: ", len(reference))

    original_seq = reference['ALIGNED_SEQ'].values[0]
    left_dna = original_seq[0:nick]
    right_dna = original_seq[nick:]
    mid_dna = ''
    check_mid = right_dna[0:4]

    if (check_mid == "TTAT"):
        right_dna = right_dna[4:]
        mid_dna = "TTAT"

    left_dna = left_dna.upper()
    mid_dna = mid_dna.upper()
    right_dna = right_dna.upper()

    rightnick =(len(left_dna) + len(mid_dna))
    return original_seq, left_dna, mid_dna, right_dna, rightnick

h_in, d_in, dir_path, nick, search_radius, min_consist_rep_size = parse_args()
original_seq, left_dna, mid_dna, right_dna, rightnick = parse_ref(h_in, nick)

# get and process data file
data_file = d_in
data_file = data_file.strip()
df = open(data_file)

base_name = os.path.basename(data_file).strip('.txt')


# get output filename
outfile = dir_path + "/" + base_name + "_consistency_log.txt"
print(outfile)
outfile = outfile.strip()
outfile = open(outfile, "w")

outfile2 = dir_path + "/" + base_name + "_consistency_table.txt"
outfile2 = open(outfile2, "w")


junction_list = []
junction_indices_list = []

linecounter = 0

for line in df:
    linecounter +=1
    line = line.strip()
    if line:
        try:
            junction_list.append(line)
            junction_indices_list.append(get_deln_indices(original_seq, line))
        except:
            print "error reading sequence in data file: line ", linecounter
            print "content of incorrect sequence is ", line
            outfile.write("error in data file on line ")
            outfile.write(str(linecounter))
            outfile.write("\n")
            outfile.write("content of incorrect sequence is ")
            outfile.write(str(line))
            outfile.write("\n")
            outfile.write("\n")

            continue

outfile.write("\n")
outfile.write("\n")


print
print "======================== RESULTS (A): ============================"
outfile.write("======================== RESULTS (A): ============================")
outfile.write("\n")
print "========= SUMMARY OF INPUT DATA AND ANALYSIS PARAMETERS =========="
outfile.write("========= SUMMARY OF INPUT DATA AND ANALYSIS PARAMETERS ==========")
outfile.write("\n")
print
outfile.write("\n")
print "%i junction sequences successfully read from file" %(len(junction_indices_list))
a = "%i junction sequences successfully read from file" %(len(junction_indices_list))
outfile.write(str(a))
print
outfile.write("\n")
z = zip(junction_list, junction_indices_list)
label = 1
table_contents = []
for pair in z:
    seq, indices = pair
    #print seq, indices
    labelstr = str(label)
    table_line = labelstr.rjust(3) + ": " + seq
    table_contents.append(table_line)
    label += 1
for line in table_contents:
    print line
    outfile.write(str(line))
    outfile.write("\n")
outfile.write("\n")
outfile.write("\n")


print
outfile.write("\n")
a = "* # bases to search left/right of junction for SD-MMEJ consistent repeats: %i" %(search_radius)
print a
outfile.write(str(a))
outfile.write("\n")
a =  "* Minimum SD-MMEJ consistent repeat size: %i" %(min_consist_rep_size)
print a
outfile.write(str(a))
outfile.write("\n")
a = "----------end summary of input data and analysis parameters--------"
print a
outfile.write(str(a))
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
print
outfile.write("\n")
a =  "======================== RESULTS (B): ============================"
print a
outfile.write(str(a))
outfile.write("\n")
a =  "=============== ALIGNMENT OF INPUT DATA SEQUENCES ================"
print a
outfile.write(str(a))
outfile.write("\n")
#get EJR objects from junction_indices_list
ejrj_list = []
for indices in junction_indices_list:
    ejrj = EJR_junction(original_seq, indices, nick, rightnick)
    ejrj_list.append(ejrj)
num_ABJs, num_MHJs = 0, 0
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ": num_ABJs +=1
    elif ejrj.jnxn_type == "MHJ": num_MHJs +=1
num_delns = len(ejrj_list)
print
outfile.write("\n")
outfile.write("\n")
a =  "* %i total junctions (%i ABJs, %i MHJs)" %(num_delns, num_ABJs, num_MHJs)
print a
outfile.write(str(a))
outfile.write("\n")
a =  "* Overhangs in original sequence in lowercase in tables below."
print a
outfile.write(str(a))
outfile.write("\n")
a =  "* Ambiguous bases at MHJs have been assigned to the left."
print a
outfile.write(str(a))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

ABJ_list = []
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ":
        ABJ_list.append(ejrj)
ori_printseq = left_dna + mid_dna.lower() + right_dna
headerstring1 = "deln".rjust(len(ori_printseq) + 7) + "  deln"
headerstring2 = "lowercase in original sequence indicates overhangs".center(len(ori_printseq)) + "bounds".center(9) + "size"
ABJ_tablestrings = []
for ejrj in ABJ_list:
    tablestring = ejrj.align_deln_raw_indices() + " " + str(ejrj.raw_deln_indices) + " " + str(ejrj.deln_size).rjust(3)
    ABJ_tablestrings.append(tablestring)
a =  "--Apparent Blunt Joins--".center(len(headerstring2))
print a
outfile.write(str(a))
outfile.write("\n")
print headerstring1
outfile.write(str(headerstring1))
outfile.write("\n")
print headerstring2
outfile.write(str(headerstring2))
outfile.write("\n")
a = "-"*len(headerstring2)
print a
outfile.write(str(a))
outfile.write("\n")
print ori_printseq
outfile.write(str(ori_printseq))
outfile.write("\n")
a  = "-"*len(headerstring2)
print a
outfile.write(str(a))
outfile.write("\n")
for tablestring in ABJ_tablestrings:
    print tablestring
    outfile.write(str(tablestring))
    outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

MHJ_list = []
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "MHJ":
        MHJ_list.append(ejrj)
MHJ_tablestrings = []
for ejrj in MHJ_list:
    tablestring = ejrj.align_deln_conv_indices() + " " + str(ejrj.converted_deln_indices) + str(ejrj.deln_size).center(6) + str(ejrj.mh_len).center(6) +str(ejrj.mh_rep_in_original_seq)
    MHJ_tablestrings.append(tablestring)
headerstring3 = " "*(len(ori_printseq) + 1) +  "deln".center(len(str(ejrj.converted_deln_indices))) + "deln".center(6) + "mh".center(6) + "mh indices".center(20)
headerstring4 = "lowercase in original sequence indicates overhangs".center(len(ori_printseq)+1) + "bounds".center(len(str(ejrj.converted_deln_indices))) + "size".center(6) + "len".center(6) + "in original seq".center(20)

a = "--Microhomology Joins--".center(len(headerstring2))
print a
outfile.write(str(a))
outfile.write("\n")

print headerstring3
outfile.write(str(headerstring3))
outfile.write("\n")

print headerstring4
outfile.write(str(headerstring4))
outfile.write("\n")

a = "-"*len(headerstring4)
print a
outfile.write(str(a))
outfile.write("\n")

print ori_printseq
outfile.write(str(ori_printseq))
outfile.write("\n")

a =  "-"*len(headerstring4)
print a
outfile.write(str(a))
outfile.write("\n")

for tablestring in MHJ_tablestrings:
    print tablestring
    outfile.write(str(tablestring))
    outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

a =  "----------end alignment of input data sequences--------"
print a
outfile.write(str(a))

print
outfile.write("\n")
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

a =  "======================== RESULTS (C): ============================"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "=============== DELETION BOUNDARY FREQUENCY DATA ================="
print a
outfile.write(str(a))
outfile.write("\n")

# maybe go back and add slice indices to tables; maybe also make tables
# that are broken out by ABJ/MHJ?
print
outfile.write("\n")
outfile.write("\n")

a =  "* Deletion boundaries are defined as the first base remaining to"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "  the left of the right-hand nick and the first base remaining to"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "  the right of the left-hand nick."
print a
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

a =  "* A consequence of this definition of deletion boundary is that"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "  values in the _del right_ and _del left_ columns in the tables below"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "  will NOT necessarily correspond to the net number of bases deleted."
print a
outfile.write(str(a))
outfile.write("\n")

a = "  See documentation for more detail as to why."
print a
outfile.write(str(a))
outfile.write("\n")

a = "* Deletion boundaries for microhomology junctions are calculated"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "  with the ambiguous bases assigned to the left of the junction."
print a
outfile.write(str(a))
outfile.write("\n")

# getting the right-hand deletion boundaries:
#put all the deln indices of the ABJs and the converted deln indices of the MHJs into a list
#the indices you want are the second ones in each tuple.
#make a list of those.
#to convert those indices to bases deleted to the right,
#find the length of the dna to the left of the left nick
#subtract that from the index, and make a new list of the
#converted values.
#make a string of the overhangs plus right dna
#for index, value in enumerate that string,
#list.count(index)
#make a new list of 3-member tuples where you glue together the index, value, and the count

indices_list = []
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ":
        indices_list.append(ejrj.raw_deln_indices)
    if ejrj.jnxn_type == "MHJ":
        indices_list.append(ejrj.converted_deln_indices)

rightbound_index_list = []
for item in indices_list:
    rightbound_index_list.append(item[1])

deln_to_right_list = []
for item in rightbound_index_list:
    bases_del_right = item - len(left_dna)
    deln_to_right_list.append(bases_del_right)

mid_and_right_dna = mid_dna + right_dna
table_line_tuple_list = []
for index, value in enumerate(mid_and_right_dna):
    freq = deln_to_right_list.count(index)
    table_line = (index, value, freq)
    table_line_tuple_list.append(table_line)

table_lines_list = []
for item in table_line_tuple_list:
    index, value, freq = item
    table_line = str(index).center(5) + str(value).center(5) + str(freq).center(5)
    table_lines_list.append(table_line)

headerstringC1 = "Deletion Boundary Frequencies"
headerstringC2 = "first base remaining to RIGHT of LEFT nick:"
headerstringC3 = "-"*15
headerstringC4= "del".center(5) + " "*10
headerstringC5 = "right".center(5) + "base".center(5) + "freq".center(5)

print headerstringC1
outfile.write(str(headerstringC1))
outfile.write("\n")

print headerstringC2
outfile.write(str(headerstringC2))
outfile.write("\n")

print headerstringC3
outfile.write(str(headerstringC3))
outfile.write("\n")

print headerstringC4
outfile.write(str(headerstringC4))
outfile.write("\n")

print headerstringC5
outfile.write(str(headerstringC5))
outfile.write("\n")

print headerstringC3
outfile.write(str(headerstringC3))
outfile.write("\n")

for line in table_lines_list:
    print line
    outfile.write(str(line))
    outfile.write("\n")
print headerstringC3
outfile.write(str(headerstringC3))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

#getting the left-hand deletion boundaries:
##put all the deln indices of the ABJs and the converted deln indices of the MHJs into a list
#the indices you want are the first ones in each tuple.
#make a list of those.
#to convert those indices to bases deleted to the left,
#find the length of left_dna + mid_dna
#subtract the index from that.
#make a new list of the # bp deleted thereby calculated
#reverse the left_dna + mid_dna string
#for index, value in enumerate that string,
#list.count(index)
#make a new list of 3-member tuples where you glue together the index, value, and the count

indices_list = []
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ":
        indices_list.append(ejrj.raw_deln_indices)
    if ejrj.jnxn_type == "MHJ":
        indices_list.append(ejrj.converted_deln_indices)

leftbound_index_list = []
for item in indices_list:
    leftbound_index_list.append(item[0])

deln_to_left_list = []
for item in leftbound_index_list:
    bases_del_left =  len(left_dna + mid_dna)  - item
    deln_to_left_list.append(bases_del_left)

left_and_center = left_dna + mid_dna
left_and_center = left_and_center[::-1]

table_line_tuple_list = []
for index, value in enumerate(left_and_center):
    freq = deln_to_left_list.count(index)
    table_line = (index, value, freq)
    table_line_tuple_list.append(table_line)

table_line_string_list = []
for item in table_line_tuple_list:
    index, value, freq = item
    table_line = str(index).center(5) + str(value).center(5) + str(freq).center(5)
    table_line_string_list.append(table_line)

headerstringC6 = "first base remaining to LEFT of RIGHT nick:"
headerstringC7 = "left".center(5) + "base".center(5) + "freq".center(5)

print headerstringC1
outfile.write(str(headerstringC1))
outfile.write("\n")

print headerstringC6
outfile.write(str(headerstringC6))
outfile.write("\n")

print headerstringC3
outfile.write(str(headerstringC3))
outfile.write("\n")

print headerstringC4
outfile.write(str(headerstringC4))
outfile.write("\n")

print headerstringC7
outfile.write(str(headerstringC7))
outfile.write("\n")

print headerstringC3
outfile.write(str(headerstringC3))
outfile.write("\n")

for line in table_line_string_list:
    print line
    outfile.write(str(line))
    outfile.write("\n")
print headerstringC3
outfile.write(str(headerstringC3))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

a =  "======================== RESULTS (D): ============================"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "============= APPARENT BLUNT JOIN SD-MMEJ ANALYSIS ==============="
print a
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

# calcs for ABJ summary
total_ABJs = 0
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ": total_ABJs +=1
num_consistent_ABJs = 0
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) != 0: num_consistent_ABJs +=1
print
outfile.write("\n")
outfile.write("\n")

a =  "----------- ABJ SUMMARY ----------"
print a
outfile.write(str(a))
outfile.write("\n")

a =  "search radius = %i" %(search_radius)
print a
outfile.write(str(a))
outfile.write("\n")

a =  "min consistent rep length = %i" %(min_consist_rep_size)
print a
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

a =  "%i total ABJs" %(total_ABJs)
print a
outfile.write(str(a))
outfile.write("\n")

a =  "%i SD-MMEJ consistent ABJs" %(num_consistent_ABJs)
print a
outfile.write(str(a))
outfile.write("\n")

a =  "----------------------------------"
print a
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

# calcs for summary table of SD-MMEJ consistent ABJs (table D1)
consistent_ABJ_list = []
consistent_ABJ_data = []
id_counter = 0
for ejrj in ejrj_list:
    id_counter += 1
    if ejrj.jnxn_type == "ABJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) != 0:
            consistent_ABJ_list.append(ejrj)
            consistent_ABJ_data.append({"id":id_counter, "type":"ABJ"})
tablestringsD1 = []
for ejrj in consistent_ABJ_list:
    tablestring = ejrj.print_jnxn().ljust(len(ori_printseq)) +str(ejrj.deln_size).center(4) + str(ejrj.raw_deln_indices).rjust(12)
    tablestringsD1.append(tablestring)
headerD1_1 = "-------------- SEQUENCES OF SD-MMEJ CONSISTENT ABJs --------------"
headerD1_2 = "Apparent point of ligation shown as upper/lowercase transition"
headerD1_3 = "-"*len(headerD1_1)
headerD1_5 = " "*len(ori_printseq) + "deln" + "deln".center(12)
headerD1_4 = "repair junction sequence".center(len(ori_printseq)) + "size" + "indices".center(12)

#print table D1
print headerD1_1
outfile.write(str(headerD1_1))
outfile.write("\n")

print headerD1_2
outfile.write(str(headerD1_2))
outfile.write("\n")

print headerD1_3
outfile.write(str(headerD1_3))
outfile.write("\n")

print headerD1_5
outfile.write(str(headerD1_5))
outfile.write("\n")

print headerD1_4
outfile.write(str(headerD1_4))
outfile.write("\n")

print headerD1_3
outfile.write(str(headerD1_3))
outfile.write("\n")

for line in tablestringsD1:
    print line
    outfile.write(str(line))
    outfile.write("\n")

print headerD1_3
outfile.write(str(headerD1_3))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

#calcs for summary table of non SD-MMEJ consistent ABJs (Table D2)
non_consistent_ABJ_list = []
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) == 0:
            non_consistent_ABJ_list.append(ejrj)
tablestringsD2 = []
for ejrj in non_consistent_ABJ_list:
    tablestring = ejrj.print_jnxn().ljust(len(ori_printseq)) +str(ejrj.deln_size).center(4) + str(ejrj.raw_deln_indices).rjust(12)
    tablestringsD2.append(tablestring)
headerD2_1 = "----------- SEQUENCES OF *NON* SD-MMEJ CONSISTENT ABJs -----------"

print headerD2_1
outfile.write(str(headerD2_1))
outfile.write("\n")

print headerD1_2
outfile.write(str(headerD1_2))
outfile.write("\n")

print headerD1_3
outfile.write(str(headerD1_3))
outfile.write("\n")

print headerD1_5
outfile.write(str(headerD1_5))
outfile.write("\n")

print headerD1_4
outfile.write(str(headerD1_4))
outfile.write("\n")

print headerD1_3
outfile.write(str(headerD1_3))
outfile.write("\n")

for line in tablestringsD2:
    print line
    outfile.write(str(line))
    outfile.write("\n")
print headerD1_3
outfile.write(str(headerD1_3))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

# ABJ aligns
for ejrj in consistent_ABJ_list:
    fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)

a =  "-------- ALIGNMENTS OF SD-MMEJ CONSISTENT ABJs ---------"
print a
outfile.write(str(a))
outfile.write("\n")
a = "--------------------------------------------------------"
print a
outfile.write(str(a))
outfile.write("\n")
j = 1
counter = 0
data_lists = []
for ejrj in consistent_ABJ_list:
    data_lists.append([])
    fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)
    t = " "*4
    a =  "** consistent ABJ %i of %i **" %(j, len(consistent_ABJ_list))
    print a
    outfile.write(str(a))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, "deletion size:", str(ejrj.deln_size)
    consistent_ABJ_data[counter]['deln'] = ejrj.deln_size
    outfile.write(str(t))
    outfile.write("deletion size:")
    outfile.write(str(ejrj.deln_size))
    outfile.write("\n")

    print t, "deletion indices: %s" %(str(ejrj.raw_deln_indices))
    outfile.write(str(t))
    outfile.write("deletion indices: %s" %(str(ejrj.raw_deln_indices)))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")

    print t, left_dna + mid_dna.lower() + right_dna, " (original sequence)"
    outfile.write(str(t))
    outfile.write(str(left_dna + mid_dna.lower() + right_dna))
    outfile.write(" (original sequence)")
    outfile.write("\n")

    print t, ejrj.align_deln_raw_indices()
    outfile.write(str(t))
    outfile.write(str(ejrj.align_deln_raw_indices()))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")

    print t, "search radius = %i" %(search_radius)
    outfile.write(str(t))
    outfile.write("search radius = %i" %(search_radius))
    outfile.write("\n")
    print t, "min consistent rep length = %i" %(min_consist_rep_size)
    outfile.write(str(t))
    outfile.write("min consistent rep length = %i" %(min_consist_rep_size))
    outfile.write("\n")
    print t, "%i total SD-MMEJ consistent repeats:" %(len(fcr_list))
    outfile.write(str(t))
    outfile.write("%i total SD-MMEJ consistent repeats:" %(len(fcr_list)))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")

    dl,dr = ejrj.raw_deln_indices
    del_len = dr - dl
    for maxpair in fcr_list:
        print t, ejrj.print_jnxn()
        outfile.write(str(t))
        outfile.write(str(ejrj.print_jnxn()))


        i1, i2, wordlen, ori = maxpair
        #adjusted wordlen,i1,i2 to account for adjacent primers
        if i1 == i2:
            a_wordlen = wordlen/2
            a_i1 = i1
            a_i2 = i2 + a_wordlen
        else:
            a_wordlen = wordlen
            a_i1 = i1
            a_i2 = i2
        data = {"len":a_wordlen}

        if a_i2 < dl:
            data["side"] = "left"
            data["rmtobr"] = rightnick - a_i1
            data["rmtodl"] = dl - a_i1
            data["p2tobr"] = rightnick - a_i2
            data["p2todl"] = dl - a_i2
            data["deltomh"] = rightnick - dr
        else:
            data["side"] = "right"
            data["rmtobr"] = a_i2 + a_wordlen - nick + del_len
            data["rmtodl"] = a_i2 + a_wordlen - dr + del_len
            data["p2tobr"] = a_i1 + a_wordlen - nick + del_len
            data["p2todl"] = a_i1 + a_wordlen - dr + del_len
            data["deltomh"] = dl - nick

        if ori == "fwd":
            data["mech"] = "loop-out"
            data["p1tobr"] = data["rmtobr"]
            data["p1todl"] = data["rmtodl"]
        else:
            data["mech"] = "snap-back"
            data["p1tobr"] = data["rmtobr"] - a_wordlen + data["p2todl"]
            data["p1todl"] = data["rmtodl"] - a_wordlen + data["p2todl"]

        if data["p1tobr"] > data["rmtobr"]:
            data["p1tobr"] = data["rmtobr"]
            data["p1todl"] = data["rmtodl"]

        data["len"] = a_wordlen
        data["p1top2"] = data["p1tobr"] - data["p2tobr"]

        motif1 = ""
        motifsnodash = ""
        for index, value in enumerate(ejrj.print_jnxn()):
            if index in range(i1, i1+wordlen):
                motif1 = motif1 + value
                motifsnodash += value
            else: motif1 = motif1 + "-"
        motif2=""
        motifsnodash += "/"
        for index, value in enumerate(ejrj.print_jnxn()):
            if index in range(i2, i2+wordlen):
                motif2 = motif2 + value
                motifsnodash += value
            else: motif2 = motif2 + "-"
        motifsnodash = motifsnodash.upper()
        if a_i2 != i2:
            motifsnodash = motifsnodash[:a_wordlen] + "/" + motifsnodash[-a_wordlen:]
        data["seq"] = motifsnodash
        data_lists[counter].append(data)
        print t, motif1, " repeated motif length:", get_motiflen(maxpair)
        outfile.write("  repeated motif length: ")
        outfile.write(str(get_motiflen(maxpair)))
        outfile.write("\n")
        outfile.write(str(t))
        outfile.write(str(motif1))
        outfile.write("\n")

        print t, motif2
        outfile.write(str(t))
        outfile.write(str(motif2))
        outfile.write("\n")

        print
        outfile.write("\n")
        outfile.write("\n")

    print headerD1_3
    outfile.write(str(headerD1_3))
    outfile.write("\n")

    j +=1
    print
    outfile.write("\n")
    outfile.write("\n")

    consistent_ABJ_data[counter]["list"] = data_lists[counter]
    counter += 1

j=1
print
outfile.write("\n")
outfile.write("\n")

outfile2.write("Sample ID\tDeletion Length\tRepair Type\tMechanism\t" +
               "Motif to Break\tMotif to Deletion\t" +
               "P1 to Break\tP1 to Deletion\tP2 to Break\t" +
               "P2 to Deletion\tP1 to P2\tMotif Length\tBreak Side\t" +
               "Deletion to MH\tMotif Sequence\n")
for elem in consistent_ABJ_data:
    for item in elem["list"]:
        outfile2.write(str(elem["id"]))
        outfile2.write("\t")
        outfile2.write(str(elem["deln"]))
        outfile2.write("\t")
        outfile2.write(elem["type"])
        outfile2.write("\t")
        outfile2.write(item["mech"])
        outfile2.write("\t")
        outfile2.write(str(item["rmtobr"]))
        outfile2.write("\t")
        outfile2.write(str(item["rmtodl"]))
        outfile2.write("\t")
        outfile2.write(str(item["p1tobr"]))
        outfile2.write("\t")
        outfile2.write(str(item["p1todl"]))
        outfile2.write("\t")
        outfile2.write(str(item["p2tobr"]))
        outfile2.write("\t")
        outfile2.write(str(item["p2todl"]))
        outfile2.write("\t")
        outfile2.write(str(item["p1top2"]))
        outfile2.write("\t")
        outfile2.write(str(item["len"]))
        outfile2.write("\t")
        outfile2.write(item["side"])
        outfile2.write("\t")
        outfile2.write(str(item["deltomh"]))
        outfile2.write("\t")
        outfile2.write(item["seq"])
        outfile2.write("\n")

a =  "------ ALIGNMENTS OF *NON* SD-MMEJ CONSISTENT ABJs --------"
print a
outfile.write(str(a))
outfile.write("\n")

a = "-----------------------------------------------------------"
print a
outfile.write(str(a))
outfile.write("\n")

for ejrj in non_consistent_ABJ_list:
    fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)
    t = " "*4
    a =  "** NON-consistent ABJ %i of %i **" %(j, len(non_consistent_ABJ_list))
    print a
    outfile.write(str(a))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, "deletion size:", str(ejrj.deln_size)
    outfile.write(str(t))
    outfile.write("deletion size:")
    outfile.write(str(ejrj.deln_size))
    outfile.write("\n")
    print t, "deletion indices: %s" %(str(ejrj.raw_deln_indices))
    outfile.write(str(t))
    outfile.write("deletion indices: %s" %(str(ejrj.raw_deln_indices)))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, left_dna + mid_dna.lower() + right_dna, "(original sequence)"
    outfile.write(str(t))
    outfile.write(str(left_dna + mid_dna.lower() + right_dna))
    outfile.write("(original sequence)")
    outfile.write("\n")
    print t, ejrj.align_deln_raw_indices()
    outfile.write(str(t))
    outfile.write(str(ejrj.align_deln_raw_indices()))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, "search radius = %i" %(search_radius)

    print t, "min consistent rep length = %i" %(min_consist_rep_size)
    outfile.write(str(t))
    outfile.write("min consistent rep length = %i" %(min_consist_rep_size))
    outfile.write("\n")

    print t, "%i total SD-MMEJ consistent repeats" %(len(fcr_list))
    outfile.write(str(t))
    outfile.write("%i total SD-MMEJ consistent repeats" %(len(fcr_list)))
    outfile.write("\n")
    a =  "---------------------------------------------------------------------"
    print a
    outfile.write(str(a))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    j+=1
a =  "----------end SD-MMEJ analysis of apparent blunt joins---------------"
print a
outfile.write(str(a))
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")




a =  "********************** RESULTS (E): ****************************"
print a
outfile.write(str(a))
outfile.write("\n")
a =  "********** MICROHOMOLOGY JUNCTION SD-MMEJ ANALYSIS *************"
print a
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

#calcs for MHJ summary
total_MHJs = 0
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "MHJ": total_MHJs +=1
num_consistent_MHJs = 0
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "MHJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) != 0:
            num_consistent_MHJs +=1

a =  "-------- MHJ SUMMARY ---------- "
print a
outfile.write(str(a))
outfile.write("\n")

a =  "search radius = %i" %(search_radius)
print a
outfile.write(str(a))
outfile.write("\n")

a =  "min consistent rep length = %i" %(min_consist_rep_size)
print a
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")

a =  "%i total MHJs" %(total_MHJs)
print a
outfile.write(str(a))
outfile.write("\n")

a =  "%i SD-MMEJ consistent MHJs" %(num_consistent_MHJs)
print a
outfile.write(str(a))
outfile.write("\n")

a =  "-------------------------------"
outfile.write(str(a))
outfile.write("\n")

print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

#calcs for summary table of SD-MMEJ consistent MHJs
consistent_MHJ_list = []
consistent_MHJ_data = []
id_counter = 0
for ejrj in ejrj_list:
    id_counter += 1
    if ejrj.jnxn_type == "MHJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) != 0:
            consistent_MHJ_list.append(ejrj)
            consistent_MHJ_data.append({"id":id_counter,"type":"MHJ"})
tablestringsE1 = []
for ejrj in consistent_MHJ_list:
    tablestring = ejrj.print_jnxn().ljust(len(ori_printseq)) +str(ejrj.deln_size).center(4) + str(ejrj.converted_deln_indices).rjust(12)
    tablestringsE1.append(tablestring)
headerE1_1 = "-------------- SEQUENCES OF SD-MMEJ CONSISTENT MHJs --------------"
headerE1_2 = "Ambiguous bases in repair product in lowercase"
headerE1_3 = "-"*len(headerE1_1)
headerE1_5 = " "*len(ori_printseq) + "deln" + "deln".center(12)
headerE1_4 = "sequence".center(len(ori_printseq)) + "size" + "indices".center(12)

#print table E1
print headerE1_1
outfile.write(str("-------------- SEQUENCES OF SD-MMEJ CONSISTENT MHJs --------------"))
outfile.write("\n")
print headerE1_2
outfile.write(str("Ambiguous bases in repair product in lowercase"))
outfile.write("\n")
print headerE1_3
outfile.write(str(headerE1_3))
outfile.write("\n")
print headerE1_5
outfile.write(str(headerE1_5))
outfile.write("\n")
print headerE1_4
outfile.write(str(headerE1_4))
outfile.write("\n")
print headerE1_3
outfile.write(str(headerE1_3))
outfile.write("\n")

for line in tablestringsE1:
    print line
    outfile.write(str(line))
    outfile.write("\n")
print headerE1_3
outfile.write(str(headerE1_3))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")


#calcs for summary table of non SD-MMEJ consistent MHJs (Table E2)
non_consistent_MHJ_list = []
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "MHJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) == 0:
            non_consistent_MHJ_list.append(ejrj)
tablestringsE2 = []
for ejrj in non_consistent_MHJ_list:
    tablestring = ejrj.print_jnxn().ljust(len(ori_printseq)) +str(ejrj.deln_size).center(4) + str(ejrj.converted_deln_indices).rjust(12)
    tablestringsE2.append(tablestring)
headerE2_1 = "----------- SEQUENCES OF *NON* SD-MMEJ CONSISTENT MHJs -----------"

print headerE2_1
outfile.write(str(headerE2_1))
outfile.write("\n")
print headerE1_2
outfile.write(str(headerE1_2))
outfile.write("\n")
print headerE1_3
outfile.write(str(headerE1_3))
outfile.write("\n")
print headerE1_5
outfile.write(str(headerE1_5))
outfile.write("\n")
print headerE1_4
outfile.write(str(headerE1_4))
outfile.write("\n")
print headerE1_3
outfile.write(str(headerE1_3))
outfile.write("\n")

for line in tablestringsE2:
    print line
    outfile.write(str(line))
    outfile.write("\n")

print headerE1_3
outfile.write(str(headerE1_3))
outfile.write("\n")


print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")

# MHJ aligns
for ejrj in consistent_MHJ_list:
    fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)

#table E3
a = "======== ALIGNMENTS OF SD-MMEJ CONSISTENT MHJs ========="
print a
outfile.write(str(a))
outfile.write("\n")
a =  "========================================================"
outfile.write(str(a))
outfile.write("\n")
j = 1
counter = 0
data_lists = []
for ejrj in consistent_MHJ_list:
    data_lists.append([])
    fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)
    t = " "*4
    a =  "** consistent MHJ %i of %i **" %(j, len(consistent_MHJ_list))
    print a
    outfile.write(str(a))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")

    print t, "deletion size: ", str(ejrj.deln_size)
    consistent_MHJ_data[counter]['deln'] = ejrj.deln_size
    outfile.write(str(t))
    outfile.write("deletion size: ")
    outfile.write(str(ejrj.deln_size))
    outfile.write("\n")

    print t, "deletion indices: %s" %(str(ejrj.converted_deln_indices))
    outfile.write(str(t))
    outfile.write("deletion indices: %s" %(str(ejrj.converted_deln_indices)))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, left_dna + mid_dna.lower() + right_dna, "  (original sequence)"
    outfile.write(str(t))
    outfile.write(str(left_dna + mid_dna.lower() + right_dna))
    outfile.write("  (original sequence)")
    outfile.write("\n")
    print t, ejrj.align_deln_conv_indices()
    outfile.write(str(t))
    outfile.write(str(ejrj.align_deln_conv_indices()))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, "search radius = %i" %(search_radius)
    outfile.write(str(t))
    outfile.write("search radius = %i" %(search_radius))
    outfile.write("\n")

    print t, "min consistent rep length = %i" %(min_consist_rep_size)
    outfile.write(str(t))
    outfile.write("min consistent rep length = %i" %(min_consist_rep_size))
    outfile.write("\n")

    print t, "%i total SD-MMEJ consistent repeats:" %(len(fcr_list))
    outfile.write(str(t))
    outfile.write("%i total SD-MMEJ consistent repeats:" %(len(fcr_list)))
    outfile.write("\n")

    print
    outfile.write("\n")
    outfile.write("\n")

    dl,dr = ejrj.converted_deln_indices
    del_len = dr - dl
    for maxpair in fcr_list:
        print t, ejrj.print_jnxn()
        outfile.write(str(t))
        outfile.write(str(ejrj.print_jnxn()))

        i1, i2, wordlen, ori = maxpair
        if i1 == i2:
            a_wordlen = wordlen/2
            a_i1 = i1
            a_i2 = i2 + a_wordlen
        else:
            a_wordlen = wordlen
            a_i1 = i1
            a_i2 = i2
        data = {"len":a_wordlen}

        if a_i2 < dl:
            data["side"] = "left"
            dl = a_i2 + a_wordlen
            data["rmtobr"] = rightnick - a_i1
            data["rmtodl"] = dl - a_i1
            data["p2tobr"] = rightnick - a_i2
            data["p2todl"] = dl - a_i2
            data["deltomh"] = rightnick - dr
        else:
            data["side"] = "right"
            dr = a_i1 + del_len
            data["rmtobr"] = a_i2 + a_wordlen - nick + del_len
            data["rmtodl"] = a_i2 + a_wordlen - dr + del_len
            data["p2tobr"] = a_i1 + a_wordlen - nick + del_len
            data["p2todl"] = a_i1 + a_wordlen - dr + del_len
            data["deltomh"] = dl - nick

        if ori == "fwd":
            data["mech"] = "loop-out"
            data["p1tobr"] = data["rmtobr"]
            data["p1todl"] = data["rmtodl"]
        else:
            data["mech"] = "snap-back"
            data["p1tobr"] = data["rmtobr"] - a_wordlen + data["p2todl"]
            data["p1todl"] = data["rmtodl"] - a_wordlen + data["p2todl"]

        if data["p1tobr"] > data["rmtobr"]:
            data["p1tobr"] = data["rmtobr"]
            data["p1todl"] = data["rmtodl"]


        data["len"] = a_wordlen
        data["p1top2"] = data["p1tobr"] - data["p2tobr"]

        motifsnodash = ""

        motif1 = ""
        for index, value in enumerate(ejrj.print_jnxn()):
            if index in range(i1, i1+wordlen):
                motif1 = motif1 + value
                motifsnodash += value
            else: motif1 = motif1 + "-"
        motif2=""
        motifsnodash += "/"
        for index, value in enumerate(ejrj.print_jnxn()):
            if index in range(i2, i2+wordlen):
                motif2 = motif2 + value
                motifsnodash += value
            else: motif2 = motif2 + "-"
        motifsnodash = motifsnodash.upper()
        if a_i2 != i2:
            motifsnodash = motifsnodash[:a_wordlen] + "/" + motifsnodash[-a_wordlen:]
        data["seq"] = motifsnodash
        data_lists[counter].append(data)

        print t, motif1, " repeated motif length: ", get_motiflen(maxpair)
        outfile.write(" repeated motif length: ")
        outfile.write(str(get_motiflen(maxpair)))
        outfile.write("\n")

        outfile.write(str(t))
        outfile.write(str(motif1))
        outfile.write("\n")

        print t, motif2
        outfile.write(str(t))
        outfile.write(str(motif2))
        outfile.write("\n")

        print
        outfile.write("\n")
        outfile.write("\n")

    print headerD1_3
    outfile.write(str(headerD1_3))
    outfile.write("\n")

    j+=1
    print
    outfile.write("\n")
    outfile.write("\n")

    consistent_MHJ_data[counter]["list"] = data_lists[counter]
    counter += 1

j=1
print
outfile.write("\n")
outfile.write("\n")

for elem in consistent_MHJ_data:
    for item in elem["list"]:
        outfile2.write(str(elem["id"]))
        outfile2.write("\t")
        outfile2.write(str(elem["deln"]))
        outfile2.write("\t")
        outfile2.write(elem["type"])
        outfile2.write("\t")
        outfile2.write(item["mech"])
        outfile2.write("\t")
        outfile2.write(str(item["rmtobr"]))
        outfile2.write("\t")
        outfile2.write(str(item["rmtodl"]))
        outfile2.write("\t")
        outfile2.write(str(item["p1tobr"]))
        outfile2.write("\t")
        outfile2.write(str(item["p1todl"]))
        outfile2.write("\t")
        outfile2.write(str(item["p2tobr"]))
        outfile2.write("\t")
        outfile2.write(str(item["p2todl"]))
        outfile2.write("\t")
        outfile2.write(str(item["p1top2"]))
        outfile2.write("\t")
        outfile2.write(str(item["len"]))
        outfile2.write("\t")
        outfile2.write(item["side"])
        outfile2.write("\t")
        outfile2.write(str(item["deltomh"]))
        outfile2.write("\t")
        outfile2.write(item["seq"])
        outfile2.write("\n")


a =  "------ ALIGNMENTS OF *NON* SD-MMEJ CONSISTENT MHJs --------"
print a
outfile.write(str(a))
outfile.write("\n")
a =  "-----------------------------------------------------------"
print a
outfile.write(str(a))
outfile.write("\n")
j=1
for ejrj in non_consistent_MHJ_list:
    fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)
    t = " "*4
    a =  "** NON-consistent MHJ %i of %i **" %(j, len(non_consistent_MHJ_list))
    print a
    outfile.write(str(a))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, "deletion size:", str(ejrj.deln_size)
    outfile.write(str(t))
    outfile.write("deletion size:")
    outfile.write(str(ejrj.deln_size))
    outfile.write("\n")
    print t, "deletion indices: %s" %(str(ejrj.raw_deln_indices))
    outfile.write(str(t))
    outfile.write("deletion indices: %s" %(str(ejrj.raw_deln_indices)))
    outfile.write("\n")
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, left_dna + mid_dna.lower() + right_dna, "(original sequence)"
    outfile.write(str(t))
    outfile.write(str(left_dna + mid_dna.lower() + right_dna))
    outfile.write("(original sequence)")
    outfile.write("\n")
    print t, ejrj.align_deln_conv_indices()
    print
    outfile.write("\n")
    outfile.write("\n")
    print t, "search radius = %i" %(search_radius)
    outfile.write(str(t))
    outfile.write("search radius = %i" %(search_radius))
    outfile.write("\n")
    print t, "min consistent rep length = %i" %(min_consist_rep_size)
    outfile.write(str(t))
    outfile.write("min consistent rep length = %i" %(min_consist_rep_size))
    outfile.write("\n")
    print t, "%i total SD-MMEJ consistent repeats" %(len(fcr_list))
    outfile.write(str(t))
    outfile.write("%i total SD-MMEJ consistent repeats" %(len(fcr_list)))
    outfile.write("\n")
    a =  "---------------------------------------------------------------------"
    print a
    outfile.write(str(a))
    outfile.write("\n")
    j+=1
    print
    outfile.write("\n")
    outfile.write("\n")
a =  "----------end SD-MMEJ analysis of microhomology joins---------------"
print a
outfile.write(str(a))
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")
print
outfile.write("\n")
outfile.write("\n")


a =  "******************************************"
print a
outfile.write(str(a))
outfile.write("\n")
a = "Analysis complete!"
print a
outfile.write(a)

print "Shown as sequence, raw deletion indices, converted deln indices, deletion size, jnxn type, #bp mh, mh in ori seq"
print "(Ambiguous bases shown in lowercase)"
print
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "MHJ":
        if len(ejrj.filtered_cr(min_consist_rep_size, search_radius)) == 0:
            consistent = "NO"
        else:
            consistent = "YES"
        print ejrj.print_jnxn(), ejrj.raw_deln_indices, ejrj.converted_deln_indices, ejrj.deln_size, "MHJ", ejrj.mh_len, ejrj.mh_rep_in_original_seq, consistent
        print

print "========= ALIGNMENTS OF CONSISTENT REPEATS: ABJs =========="
print
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "ABJ":
        fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)
        print "============================================================"
        print "For deletion indices " + str(ejrj.raw_deln_indices) + ejrj.jnxn_type + ":"
        print left_dna + mid_dna.lower() + right_dna, "(original sequence)"
        print ejrj.align_deln_raw_indices(), ejrj.raw_deln_indices
        print
        print "With minimum consistent rep length %i, search radius %i:" %(min_consist_rep_size, search_radius)
        print "%i total SD-MMEJ consistent repeats: %s" %(len(fcr_list), str(fcr_list))
        print
        for maxpair in fcr_list:
            print ejrj.print_jnxn()
            i1, i2, wordlen, ori = maxpair
            motif1 = ""
            for index, value in enumerate(ejrj.print_jnxn()):
                if index in range(i1, i1+wordlen):
                    motif1 = motif1 + value
                else: motif1 = motif1 + "-"
            motif2=""
            for index, value in enumerate(ejrj.print_jnxn()):
                if index in range(i2, i2+wordlen):
                    motif2 = motif2 + value
                else: motif2 = motif2 + "-"
            print motif1, maxpair, "motif length is", get_motiflen(maxpair)
            print motif2
            print
    print

print "========= ALIGNMENTS OF CONSISTENT REPEATS: MHJs =========="
print """
    ***NOTE***: for purposes of identifying individual deletion events,
    deletions are given as generated by the code that determines the set
    of valid deletion events.  Ambiguous bases are assigned to the left
    in the alignments with the original sequences, and the resulting
    converted deletion boundaries indicated next to the alignment.
    Note that this results in the some initially distinct
    deletion events being indistinguishable in the alignment.
    """
print
for ejrj in ejrj_list:
    if ejrj.jnxn_type == "MHJ":
        fcr_list = ejrj.filtered_cr(min_consist_rep_size, search_radius)
        print "============================================================"
        print "For deletion indices " + str(ejrj.raw_deln_indices) + ejrj.jnxn_type + ":"
        print left_dna + mid_dna.lower() + right_dna, "(original sequence)"
        print ejrj.align_deln_conv_indices(), ejrj.converted_deln_indices
        print
        print "With minimum consistent rep length %i, search radius %i:" %(min_consist_rep_size, search_radius)
        print "%i total SD-MMEJ consistent repeats: %s" %(len(fcr_list), str(fcr_list))
        print
        for maxpair in fcr_list:
            print ejrj.print_jnxn()
            i1, i2, wordlen, ori = maxpair
            motif1 = ""
            for index, value in enumerate(ejrj.print_jnxn()):
                if index in range(i1, i1+wordlen):
                    motif1 = motif1 + value
                else: motif1 = motif1 + "-"
            motif2=""
            for index, value in enumerate(ejrj.print_jnxn()):
                if index in range(i2, i2+wordlen):
                    motif2 = motif2 + value
                else: motif2 = motif2 + "-"
            print motif1, maxpair, "motif length is", get_motiflen(maxpair)
            print motif2
            print
print

print "exiting SD-MMEJ data analysis program..."
outfile.close()
outfile2.close()















