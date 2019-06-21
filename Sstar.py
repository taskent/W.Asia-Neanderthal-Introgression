##S* statistics 
##author: Recep Ozgur Taskent
##date: 06.19.2017
##run as python Sstar.py inputfile outputfile outputfile2

import itertools
from itertools import islice
import random
import sys
import os

inputfile = sys.argv[1]  #input file showing genotypes for each individual (in .bed format)
outputfile = sys.argv[2]  #output file for the haplotypes detected by S* 
outputfile2 = sys.argv[3] #output file for the total number of segregating sites used for S* calculations



##Calculate 50kb sliding windows with 20kb interval 
inf = (inputfile, 'r')
len_inf = sum(1 for line in inf) #length of the input file 
inf = open(inputfile, 'r')
snp_pos_list = [] 
count = 0
for line in inf: #get the positions of the first and the last SNP in the data set   
    count += 1
    l = line.split()
    chr_no  = l[0]
    snp_pos = int(l[2])
    snp_pos_list.append(snp_pos)
    if count == 1:
        start_snp = snp_pos
    elif count == len_inf:
        stop_snp = snp_pos

snp_len = int(stop_snp - start_snp) #total size of the

numOfChunks = ((snp_len-50000)/20000)+1   #number of 50kb chunks/windows with 20kb sliding interval
start_pos_list = []
stop_pos_list = []
for i in range(numOfChunks):  #save the position of the first SNP for each window 
    start_pos = start_snp + (int(i)*20000)
    start_pos_list.append(start_pos)
for i in range(numOfChunks): #save the position of the last SNP for each window
    stop_pos = start_snp + 50000 + (int(i)*20000)
    stop_pos_list.append(stop_pos)
    




wind_start_snp_dict = {}   #create an empty dictionary for the starting rows in the input file for each window  
wind_stop_snp_dict = {}	   #create an empty dictionary for the ending rows in the input file for each window


for j in range(numOfChunks):   #for each window
    inf = open(inputfile, 'r')
    j = int(j)
    mini = float('inf')   #start a variable for the starting row number (the lowest row number) for each window 
    maxi = 0			  #start a variable for the ending row number (the highest row number) for each window
    count6 = 0			  #start a counter for the row number in the input file
    for line in inf:  
        count6 += 1  	#increase the row number by one
        line = line.split()
        pos = int(line[2])
        if (pos >= int(start_pos_list[j])) and (pos <= int(stop_pos_list[j])):  #if the SNP position is in between the starting and ending positions of that window:
            if count6 > maxi:   	#if the row number counter is bigger than maxi (variable for the highest row number in the window), assign maxi to the count 
                maxi = count6
            if count6 < mini: 		#if the row number counter is smaller than mini (variable for the lowest row number in the window), assign mini to the count 
                mini = count6
    if mini != float('inf') and maxi != 0:   #if for the windows there are more than one SNP, that is, if there are some SNPs within the given window, add the lowest and highest row numbers for that window to the corresponding dictionaries  
        wind_start_snp_dict[j] = mini
        wind_stop_snp_dict[j] = maxi
    





##Indices of samples belonging to different populations in 1000 Genomes Phase I data set 
gbr = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,358]
fin = [56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181]
chs = [182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287]
pur = [226,227,260,261,262,263,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347]
clm = [322,323,324,325,326,327,328,329,330,331,332,348,349,350,351,352,353,354,355,356,357,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,402,403]
ibs = [396,397,398,399,400,401,404,405,406,407,408,409,410,411]
ceu = [412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496]
Yri = [497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778]
chb = [516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612]
jpt = [634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730]
lwk = [693,694,695,696,697,698,699,700,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,859,860,861,862,863,864,865,866,867]
asw = [868,888,889,890,891,892,893,894,895,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994]
mxl = [869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942]
tsi = [995,996,997,998,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091,1092]




for i in chb:  			
    i = int(i) + 4
    e_asn.append(i) 	#Columns for individuals with East Asian individuals in 1000G phase I data set - in .bed format
for i in chs:			
    i = int(i) + 4
    e_asn.append(i)		#Columns for individuals with East Asian ancestry in 1000G phase I data set - in .bed format
for i in jpt:			
    i = int(i) + 4
    e_asn.append(i)		#Columns for individuals with East Asian ancestry in 1000G phase I data set - in .bed format
for i in gbr:			
    i = int(i) + 4
    nw_eur.append(i)	#Columns for individuals with Northern and Western European ancestry in 1000G phase I data set - in .bed format
for i in fin:
    i = int(i) + 4
    nw_eur.append(i)	#Columns for individuals with Northern and Western European ancestry in 1000G phase I data set - in .bed format
for i in ceu:
    i = int(i) + 4
    nw_eur.append(i)	#Columns for individuals with Northern and Western European ancestry in 1000G phase I data set - in .bed format
yri = []
for i in Yri:
    i = int(i) + 4
    yri.append(i)		#Columns for individuals from Yoruba in 1000G phase I data set - in .bed format






outf = open(outputfile, 'w')  #open an output file for the haplotypes detected by S*
outf2 = open(outputfile2, 'w')  #open an output file for the total number of segregating sites used for S* calculations

snp_geno_dict = {} 	#set an empty snp-genotype dictionary
snp_pos_dict = {}	#set an empty snp-position dictionary
snp_list = [] 	#set an empty snp list	
snp_list_wind_sstats = []	#set an empty list for derived snps where Sstats can be calculated - equivalent of mutation rate in the ms simulations
haplo_dict = {}		#set an empty dictionary for haplotypes
snp_count = 0	 #set a counter for the snps where S*stats will be calculated
count_totalsnp_wind = 0		#set a counter for snps where 20 Northern and Western Europe individuals have at least one derived allele  
count_totalsnp_sstats = 0	#set a counter for snps for which S*stats can be calculated - equivalent of mutation rate in ms simulations
haplo_list = []		#set an empty list for haplotypes

for j in range(len(wind_start_snp_dict)):  ##for each window with the ascending order:
	j = int(j) 
	count_j = 0
	for KEY in wind_start_snp_dict.keys():
			if KEY == j:
					count_j = 1
	if count_j == 1:
			pass
	elif count_j == 0:
			continue
	yri13 = random.sample(yri, 13) #get random 13 individuals from Yoruba
	nw_eur20 = random.sample(nw_eur, 20)	#get random 20 individuals from Northern and Western Europe
	for m in nw_eur20: 	#for all 20 individuals from Northern and Western Europe
		inF = open(inputfile, 'r') 		#Read the inputfile  
		prev_snp = 0 	#set the position for the previous snp to zero, to change it later
		for LINE in islice(inF, int(wind_start_snp_dict[j]) - 1, int(wind_stop_snp_dict[j])):: 	#for each line in the input file within given window 
			l = LINE.split()	
			if int(l[m]) > 0: 	#if there is at least one derived allele at the current snp position for the individual for which S*score is going to be calculated, otherwise, skip the position
				count_totalsnp_wind += 1	#increase the total number of snps for which there is a derived allele for the 20 individuals from Yoruba
				check = 0 		#start a check vector for 13 Yoruba individuals - to check whether they have a derived allele at the current snp position, if so, skip the position
				if int(l[2]) - prev_snp >= 10 and (l[2] not in snp_list_wind_sstats):		#if the corrent snp position is not within 10 bp of the previous snp position and the corrent snp position is not in the list of snps for which S*stats have been calculated for other individuals
					count_totalsnp_sstats += 1		#increase the total number of derived snp in the window for which S*stats can be calculated - equivalent of mutation rate in ms simulations
					snp_list_wind_sstats.append(l[2])	#add current snp position to the list for derived snps where Sstats can be calculated - equivalent of mutation rate in ms simulations
					pass
				else:	
					continue
				for k in yri13:		#check whether any of 13 Yoruba individuals has a derived allele at the current snp position, if so, skip the position
					if l[k] != "0":	
						check = 1
				if check == 1:	
					continue
				snp_count += 1		#increase snp count by one
				snp_list.append(snp_count)	 #add snp count to the list for the snps for which S*stats will be calculated   
				geno_list = []		#create an empty genotype list for 20 Northern and Western European individuals for the current snp
				for k in nw_eur20:
					geno_list.append(l[k])		#add genotypes of 20 Northern and Western European individuals for the current snp to the genotype list
				snp_geno_dict[snp_count] = geno_list	#add genotype list to the snp-genotype dictionary for the current snp 
				snp_pos_dict[snp_count] = int(l[2])		#add current snp position to the snp-position dictionary for the current snp 
				prev_snp = int(l[2])	#reset previous snp position to the current snp position
			elif int(l[m]) == 0:	#if there is no derived snp at the current snp position for the Yoruba individual, but if any of the 13 Northern and Western Europeans has at least one derived allele...
				for k in nw_eur13:	
					if l[k] != "0" and (l[2] not in snp_list_wind_sstats):
						count_totalsnp_sstats += 1		#increase the total number of derived snp in the window for which S*stats can be calculated - equivalent of mutation rate in ms simulations  
						snp_list_wind_sstats.append(l[2])		#add current snp position to the list for derived snps where S*stats can be calculated - equivalent of mutation rate in ms simulations                	
		inF.close()		#close the file
		score_G_x_dict = {} #empty dictionary for snp scores
		score_x_y_dict = {}	#empty dictionary for snp-pair scores
		score_G_x = 0
		for x in range(len(snp_list)-1):   #for each snp and snp pair, start S*score vectors from zero 
			for y in range(x+1, len(snp_list)):
				score_G_x_dict[snp_list[x]] = 0
				score_G_x_dict[snp_list[y]] = 0
		for x in range(len(snp_list)-1):   #from now on, calculate S*stats
			for y in range(x+1, len(snp_list)):
				total_mis = 0 	#start a mismatch vector for the snp-pair
				e = snp_list[x] 	#first snp position in the file - particular row where the snp is located in the file
				f = snp_list[y] 	#second snp position in the file - particular row where the snp is located in the file
				geno_list_snp1 = snp_geno_dict[e]  #genotypes for all 20 Northern and Western Europeans for the first snp
				geno_list_snp2 = snp_geno_dict[f]  	#genotypes for all 20 Northern and Western Europeans for the second snp
				for v in range(len(geno_list)): 	#for each individual from nweur get the genotypes  
					w = geno_list_snp1[v] 
					z = geno_list_snp2[v]
					total_mis += abs(int(w) - int(z)) 	#for the snp-pair, calculate the total mismatch between the genotypes of each individual from 20 total individuals from yoruba  
				if total_mis == 0: 		#if there is no mismatch for the genotypes of all 20 individuals for the two snps      
					if score_G_x_dict[e] != -1 * float('inf'):  	#if S(G, x) is not minus infinity, that is, if the first snp (e) does not have a minus infinity score retained from earlier snp-pair comparisons 
						S_x_y = 5000 + int(snp_pos_dict[f])- int(snp_pos_dict[e]) 	#instantaneous score for the snp-pair is: 5000 + the nucleotide distance between the two snps
						score = S_x_y + score_G_x_dict[e] 	#cumulative S*score for the snp pair is:  instantaneous score for the snp pair + the first snp's retaining individual score
						if score_G_x_dict[f] < score: 	#if the retaining individual score for the second snp is smaller than the cumulative score of the snp pair, equate it to score, that is, if the second snp has a better score for this snp-pair comparison, equate its individual retaining score to the cumulative score of the snp pair   
							score_G_x_dict[f] = score 
						snp_pair = str(e) + "_" + str(f) 
						score_x_y_dict[snp_pair] = score  	#add cumulative score of the snp pair to the snp-pair score dictionary
					elif score_G_x_dict[e] == -1 * float('inf'):   #if S(G, x) is minus infinity, that is, if the first snp (e) does have a minus infinity score retained from earlier snp pair   
						S_x_y = 5000 + int(snp_pos_dict[f])- int(snp_pos_dict[e])	#instantaneous score for the snp pair is: 5000 + the nucleotide distance between the two snps
						score = S_x_y	#cumulative S*score for the snp-pair is:  instantaneous score for the snp pair 
						if score_G_x_dict[f] < score: 	#if the retaining individual score for the second snp is smaller than the cumulative score of the snp pair, equate it to score, that is, if the second snp has a better score for this snp-pair comparison, equate its individual retaining score to the cumulative score of the snp pair
							score_G_x_dict[f] = score
						snp_pair = str(e) + "_" + str(f)
						score_x_y_dict[snp_pair] = score 	#add cumulative score of the snp pair to the snp-pair score dictionary
				elif total_mis > 0 and total_mis < 6: #if there are 1 to 5 mismatches btw two snps
					if score_G_x_dict[e] != -1 * float('inf'):   #if S(G, x) is not -infinity, that is, if the first snp (e) does not have a minus infinity score retained from earlier snp pair comparisons 
						S_x_y = int(-1 * 10000) + int(snp_pos_dict[f]) - int(snp_pos_dict[e])	#instantaneous score for the snp pair is: -10000 + the nucleotide distance between the two snps
						score = S_x_y + score_G_x_dict[e]	#cumulative S*score for the snp pair is:  instantaneous score for the snp pair + the first snp's individual retaining score
						if score_G_x_dict[f] < score:	#if the retaining individual score for the second snp is smaller than the cumulative score of the snp pair, equate it to score, that is, if the second snp has a better score for this snp-pair comparison, equate its individual retaining score to the cumulative score of the snp pair
							score_G_x_dict[f] = score
						snp_pair = str(e) + "_" + str(f)
						score_x_y_dict[snp_pair] = score	#add cumulative score of the snp pair to the snp-pair score dictionary
					elif score_G_x_dict[e] == -1 * float('inf'):   #if S(G, x) is -infinity, that is, if the first snp (e) does have a minus infinity score retained from earlier snp-pair comparisons 
						S_x_y = int(-1 * 10000) + int(snp_pos_dict[f]) - int(snp_pos_dict[e])	#instantaneous score for the snp pair is: -10000 + the nucleotide distance between the two snps
						score = S_x_y	#cumulative S*score for the snp pair is:  instantaneous score for the snp-pair 
						if score > 0:	#if score is positive, equate it to the retaining individual score for the second snp 
							score_G_x_dict[f] = score 
						snp_pair = str(e) + "_" + str(f)
						score_x_y_dict[snp_pair] = score	#add cumulative score of the snp pair to the snp-pair score dictionary
				elif total_mis > 5:  #if there are more than 5 mismatches between the two snps
					score = -1 * float('inf') 	#equate the cumulative score for the snp pair to minus infinity
					snp_pair = str(e) + "_" + str(f)
					score_x_y_dict[snp_pair] = float(score) #add cumulative score of the snp-pair to the snp-pair score dictionary
		maxi = 0
		for snp in score_G_x_dict.keys(): 	#get snp with the maximum S*score  
			if score_G_x_dict[snp] > maxi:
				maxi = float(score_G_x_dict[snp])
				maxScore_snp = str(snp)
		maxi = 0
		for snp_pair in score_x_y_dict.keys(): 	 #get snp-pair with the maximum S*score
			if score_x_y_dict[snp_pair] > maxi:
				maxi = float(score_x_y_dict[snp_pair])
				maximumScore = float(score_x_y_dict[snp_pair])
				maxScore_snp_pair = str(snp_pair)
		if maxi > 0: 	#if there is a snp-pair with >0 S*score 
			haplo_list = []		##add all snps that contributed to the maximum S*score to the haplotype list, start with the maximum scored snp pair...
			x_snp = maxScore_snp_pair.split('_')[0]		#first snp in the maximum scored snp pair
			y_snp = maxScore_snp_pair.split('_')[1]		#second snp in the maximum scored snp pair
			y_snp_pos = snp_pos_dict[int(y_snp)] 	#get second snp's position
			haplo_list.append(y_snp_pos)  ##add the position for the second snp to the haplotype list 
			while True:  ##add all the snps within the haplotype, continue with the other snps within the haplotype
				maximum = 0
				check = 0
				count5 = 0
				if x_snp == "1":  ##if the first snp in the max scored snp pair is the first snp within the window, add that snp and break...
					x_snp_pos = snp_pos_dict[int(x_snp)]
					haplo_list.append(x_snp_pos)
					break
				for snp_pair in score_x_y_dict.keys(): ##otherwise, for each snp pair for which S*score was calculated, look for the first snp in the maximum scored snp pair...
					first_snp = snp_pair.split("_")[0]		#first snp in the snp pair
					second_snp = snp_pair.split("_")[1]	#second snp in the snp pair
					if second_snp == x_snp: 	#if second snp in the current snp pair is the first snp in the maximum scored snp pair  
						if score_x_y_dict[snp_pair] > maximum:  #and if the S*score of the current snp pair > 0
							maximum = score_x_y_dict[snp_pair] 	#set maximum to S*score of the current snp pair 
							maxScore_snp_pair = str(snp_pair) 	#set maximum-scored snp pair to the current snp pair to get retrospectively all snps that are linked to the first maximum-scored snp pair
							count5 = 1		#set count vector to 1
				if count5 == 1 and maximum > 0: 	#if count is equal to one, that is, if there are more than snp-pair included in the haplotype
					x_snp = maxScore_snp_pair.split('_')[0]		#first snp in the maximum scored snp pair
					y_snp = maxScore_snp_pair.split('_')[1]		#second snp in the maximum scored snp pair
					y_snp_pos = snp_pos_dict[int(y_snp)]		#get second snp's position
					haplo_list.append(y_snp_pos)		##add the position for the second snp to the haplotype list
				if maximum == 0:	#if maximum is equal to zero, that is, if there is no other snp pair in the haplotype than the first maximum-scored snp pair  
					x_snp_pos = snp_pos_dict[int(x_snp)]	#get first snp's position
					haplo_list.append(x_snp_pos)	##add the position for the first snp (from the first maximum scored snp-pair) to the haplotype list and break
					break
			haplo_dict[j] = haplo_list
		snp_geno_dict = {} #directory for genotypes at each snp. Save as `snp_count` : `genotype_list`.
		snp_list = [] #list for snps (the number of snp where the individual of concern is polymorphic while the reference individual(s) carry the ancestral alleles).
		snp_pos_dict = {} #dictionary for snp positions. save as `snp_count` : `snp_position`.
		count4 = 0
		if not haplo_list == []:
			highest = int(haplo_list[0]) 	#the end position of the haplotype 
			lowest = int(haplo_list[len(haplo_list)-1]) #the start position of the haplotype
			length = highest -lowest  #length of the haplotype
			outf.write("%s\t%s\t%f\t%s\t%d\n" % (j, m, maximumScore, haplo_list, length)) #write the window number, individual for which the S* was calculated, S* of the haplotype to the output file,  the haplotype, and the length of the haplotype to the outputfile
		haplo_list = []
	outf2.write("%s\t%d\n" % (j, count_totalsnp_sstats))  #write the window number and the total number of segregating sites used for S* calculations to the outputfile2
	snp_list_wind_sstats = []
	count_totalsnp_sstats = 0
	count_totalsnp_wind = 0	
outf.close()    #close outputfile
outf2.close()	#close outputfile2
		
		
