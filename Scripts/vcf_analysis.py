import sys
from cyvcf2 import VCF
from collections import defaultdict
import collections
import pandas as pd



def main():
    vcf_file=sys.argv[1]
    outfile=sys.argv[2]
    outfile_del=sys.argv[3]
    len_dict = defaultdict(list)
    row=[]
    complete=[]
    mate_id=[]
    del_row=[]
    del_mat=[]
    for variant in VCF(vcf_file):
        if(variant.FILTER!='MinQUAL' and variant.FILTER!='MinGQ' and variant.INFO.get('PRECISE')):
            #print(variant.FILTER)
            #break
            #print(variant.INFO.get('SVTYPE'),variant.CHROM,variant.INFO.get('MATEID'),variant.ID,variant.POS)
            if(variant.INFO.get('SVTYPE') == 'DEL'):
            #print(variant)
                #if(int(variant.INFO.get('SVLEN'))<-1000):
                #print(int(variant.INFO.get('SVLEN')))
                #break
                        
                del_row.append(variant.INFO.get('SVTYPE'))
                del_row.append(variant.CHROM)
                del_row.append(variant.INFO.get('SVLEN'))
                del_row.append(variant.POS)
                del_row.append(variant.REF)
                del_row.append(variant.ALT)
            #break
            if (variant.INFO.get('SVTYPE') == 'BND'):
                print (variant)
                #print(variant.INFO.get('SVTYPE'),variant.CHROM,variant.INFO.get('MATEID'),variant.ID,variant.POS)
                #break
                #len_dict[variant.INFO.get('SVTYPE')].append(variant.INFO.get('MATEID'))
            
                row.append(variant.INFO.get('SVTYPE'))
                row.append(variant.CHROM)
                row.append(variant.POS)
                row.append(variant.REF)
                row.append(variant.ALT)
                row.append(variant.INFO.get('RE'))
    i=0
    while i<len(row):
        complete.append(row[i:i+6])
        i+=6

    df = pd.DataFrame(complete)
    
    df.columns=['SVTYPE','chr','POS','REF','ALT','No.Of.Reads']
    
    print(df.size)  
    #translocations=extract_translocations(df)
    df.to_csv(outfile)
    
    
    #print(translocations)
    #df1=df[df['ID'].isin(translocations)]
    #print(df)
    #for index, row in df.iterrows():
        #print(df['ID'])
        #if(df['ID'] in (ID)):
            #print(row) 
    #print(MATEID)
    #df1=df.isin({'mateid':ID}) 
    #df1.to_csv('BND1.csv') 
    #print(df1)
    #print(df.loc[df1['mateid'] == 'True'])
    #print(df.isin({'mateid':MATEID}))
    i=0
    while i<len(del_row):
        del_mat.append(del_row[i:i+6])
        i+=6

    df_del = pd.DataFrame(del_mat)
    df_del.columns=['SVTYPE','chr','Len','POS','REF','ALT']
    df_del.to_csv(outfile_del)
    
def extract_translocations(df):
    trans = pd.DataFrame(columns=['SVTYPE','chr','mateid','ID','POS','REF','ALT','SVTYPE_trans','chr_trans','mateid_trans','POS_trans','REF_trans','ALT_trans'])
    print(df)
    ID= list(df['ID'])
    #print(len(ID))
    MATEID=list(df['mateid'])
    for element in ID:
        row1=df.loc[df['ID'] == element] 
        row2=df.loc[df['mateid'] == element]
        row2.columns=['SVTYPE_trans','chr_trans','ID','mateid_trans','POS_trans','REF_trans','ALT_trans']
        #print(row1,row2)
        if (str(row1['chr'])!=str(row2['chr_trans'])):
            #print('true')
            row1=row1.join(row2.set_index('ID'), on='ID')
            #print(row1)
            #break
            trans=trans.append(row1,ignore_index=True)
            #trans=trans.append(row2,ignore_index=True)
            df=df[df['ID'] != element]
            df=df[df['mateid'] != element]
            #print(trans)
    return trans
            

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 


if __name__ == '__main__':
    main()
            



