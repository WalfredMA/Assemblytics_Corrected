#! /usr/bin/env python
import argparse

import gzip
# Author: Maria Nattestad
# Email: mnattest@cshl.edu
# This script is part of Assemblytics, a program to detect and analyze structural variants from an assembly aligned to a reference genome using MUMmer. 
    
    


def run(args):
    filename = args.delta
    minimum_variant_size = args.minimum_variant_size
    #minimum_variant_size=1
    f = open(filename)
    header1 = f.readline()
    if header1[0:4]=="\x1f\x8b\x08\x08":
        f.close()
        f = gzip.open(filename)
        header1 = f.readline()
    
    # Ignore the first two lines for now
    f.readline()

    linecounter = 0

    current_reference_name = ""
    current_reference_position = 0

    current_query_name = ""
    current_query_position = 0
    
    std='+'
    std0=1

    variants = []
    
    distances = []
    
    combine_num=0

    for line in f:
        if line[0]==">":
            # linecounter += 1
            # if linecounter > 1:
            #     break
            fields = line.strip().split()
            current_reference_name = fields[0][1:]
            current_query_name = fields[1]
            
        else:
            fields = line.strip().split()
            if len(fields) > 4:
                
                if int(fields[2]) <= int(fields[3]):
                
                    std='+'
                    std0=1
                    # current_reference_position = int(fields[0])
                    current_reference_position = min(int(fields[0]),int(fields[1]))
                    # fields[1] is the reference position at the end of the alignment
                    # current_query_position = int(fields[2])
                    current_query_position = min(int(fields[2]),int(fields[3]))
                    # fields[3] is the query position at the end of the alignment
                else:
                    
                    std='-'
                    std0=-1
                    current_reference_position = min(int(fields[0]),int(fields[1]))
                    current_query_position = max(int(fields[2]),int(fields[3]))
            
                report = [current_reference_name,current_reference_position,current_reference_position,"Assemblytics_w_"+str(len(variants)+1),0,std,"",0,0,current_query_name,"within_alignment",current_query_position,current_query_position]        
                distances=[] 
      
                mis=0
                match=0

            else:

                tick = int(fields[0])
                size = abs(tick)
                allow_distance=2*(report[7]+report[8])-report[4]

                if size>0 and size <= min(allow_distance+2,19): # then go back and edit the last entry to add 1 more to its size
                    #report = variants[-1]
                    if size!=1:
                        combine_num+=1
                        
                    mis=mis+1
                    match=match+size-1
                    report[4] = report[4] + size # size
                    if tick > 0: # deletion, moves in reference
                        report[7] = report[7] + 1  # reference gap size
                        
                        report[2] = report[2] + size # reference end position
                        current_reference_position += size # update reference position after deletion
                        
                        report[12] = report[12] + std0*(size-1)
                        current_query_position += std0*(size-1)
                    elif tick < 0: # insertion, moves in query
                        report[8] = report[8] + 1 # query gap size
                        
                        report[2] = report[2] + size-1
                        current_reference_position += size-1 # query end position
                        
                        report[12] = report[12] + std0*size
                        current_query_position += std0*size # update query position after insertion
                else: # report the last one and continue


                    newdistances=[x for x in distances]
                    
                    for distance,score_record in distances:

                        allow_distance0=2*(report[7]+report[8])-report[4]+score_record

                        if  distance>min(allow_distance0,19):

                            break

                        else:
                            
                            combine_num+=1

                            combined_report=variants[-1]

                            report[1]=combined_report[1]

                            report[3]=combined_report[3]
                            
                            report[4]=combined_report[4]+report[4]+distance

                            report[7]=combined_report[7]+report[7]

                            report[8]=combined_report[8]+report[8]			

                            report[11]=combined_report[11]

                            variants=variants[:-1]

                            newdistances=newdistances[1:]							

                    distance=newdistances

                    if report[4]-max(report[7],report[8])>5:
                        
                        if report[7]>report[8]:
                            
                            report[6]="Repeat_contraction"
                        
                        else:
                            
                            report[6]="Repeat_expansion"
                            
                    else:
                        
                        if report[7]>report[8]:
                            
                            report[6]="Deletion"
                        
                        else:
                            
                            report[6]="Insertion"
                        
                    
                    
                    variants.append(report)


                    distances=[(abs(tick)-1,2*(report[7]+report[8])-report[4])]+distances

                
                    mis=0
                    match=0
                    
                    current_reference_position += abs(tick) -1
                    current_query_position += std0*(abs(tick) -1) 
                    if tick > 0:
                        
                        # report = "%s\t%d\t%d\tAssemblytics_%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (current_reference_name,current_reference_position,current_reference_position+size,len(variants)+1,size,"+","Deletion",size,0,current_query_name,"within_alignment")
                        #report = [current_reference_name,current_reference_position,current_reference_position+size,"Assemblytics_w_"+str(len(variants)+1),size,std,"Deletion",size,0,current_query_name,"within_alignment",current_query_position,current_query_position, current_query_position+std0*size]
                        
                        report = [current_reference_name,current_reference_position,current_reference_position+1,"Assemblytics_w_"+str(len(variants)+1),1,std,"",1,0,current_query_name,"within_alignment",current_query_position,current_query_position] 
                        
                        current_reference_position += 1 # update reference position after deletion
                        
                    elif tick < 0:
                        
                        
                        # report = "%s\t%d\t%d\tAssemblytics_%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (current_reference_name,current_reference_position,current_reference_position,len(variants)+1,size,"+","Insertion",0,size,current_query_name,"within_alignment")
                        #report = [current_reference_name,current_reference_position,current_reference_position,"Assemblytics_w_"+str(len(variants)+1),size,std,"Insertion",0,size,current_query_name,"within_alignment",current_query_position,current_query_position+size]
                         # update query position after insertion
                   
                        report = [current_reference_name,current_reference_position,current_reference_position,"Assemblytics_w_"+str(len(variants)+1),1,std,"",0,1,current_query_name,"within_alignment",current_query_position,current_query_position+std0*1] 
                        
                        current_query_position += std0*1 # update query position after insertion
                        
                # TESTING
                # print line, report
                

    f.close()
    

    newcounter = 1
    for line in variants:
        # report = "%s\t%d\t%d\tAssemblytics_%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % line
        if line[4] >= minimum_variant_size:
            line[3] = "Assemblytics_w_%d" % (newcounter)
            
            print "\t".join(map(str,line[0:10])) + ":" + str(line[11]) + "-" + str(line[12]) + ":+"+"\t" + line[10]
            # print "\t".join(map(str,line))
            newcounter += 1


def main():

    parser=argparse.ArgumentParser(description="Outputs MUMmer coordinates annotated with length of unique sequence for each alignment")
    parser.add_argument("--delta",help="delta file" ,dest="delta", type=str, required=True)
    parser.add_argument("--min",help="Minimum size (bp) of variant to include, default = 50" ,dest="minimum_variant_size",type=int, default=50)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()




