###############################################################
# dNdS statistic simple model (equal nucleotide change raiting)
# parameter: name=M (ancestor strain)
#            name_ref=N (reference genome name)
###############################################################
{
	if(FNR==1)
		num_file++;
	if(num_file==1)
	{
		arr_prob[$3,$4,1]=$5;
		arr_prob[$3,$4,2]=$6;
		arr_prob[$3,$4,3]=$7;
		arr_prob[$3,$4,4]=$8;
		arr_prob[$3,$4,5]=$9;
		if($7==1)
		{
			arr_prob_1[$3,1]=$5;
			arr_prob_1[$3,2]=$6;
			arr_prob_1[$3,3]=$7;
		}
	}
	if(num_file==2)
	{
		if(gene!=$2)
		{
			if(FNR>1)
			{
				len=length(arr_seq[1])
				for(start=1;start<=len&&substr(arr_seq[num_anc],start,1)=="-";start++)
					num_nuc+=(substr(arr_seq[num_ref],start,1)!="-") ? 1 : 0;
				for(;start<=len&&num_nuc%3!=0;start++)
					num_nuc+=(substr(arr_seq[num_ref],start,1)!="-") ? 1 : 0;
#				print substr(arr_seq[num_anc],start);
#				print substr(arr_seq[num_ref],start);

				for(i=1;i<=pos;i++)
				{
					num_N=num_S=sum_N=sum_S=num=num_all=sum_all=0;
					if(i!=num_anc)
					{
						seq_anc=seq_comp=seq_ref="";
						for(j=start;j<=len;j++)
						{
							if(substr(arr_seq[num_ref],j,1)!="-")
							{
								seq_anc=seq_anc substr(arr_seq[num_anc],j,1);
								seq_comp=seq_comp substr(arr_seq[i],j,1);
								seq_ref=seq_ref substr(arr_seq[num_ref],j,1);
							}
						}
#						print gene, num_anc, seq_anc;
#						print gene, i, seq_comp;
#						print gene, num_ref, seq_ref;
						len_seq=length(seq_ref);
						
						for(j=1;j<=len_seq;j+=3)
						{
							num_all++;
							cod_anc=substr(seq_anc,j,3);
							cod_prob=substr(seq_comp,j,3);
							if(cod_anc!=cod_prob&&arr_prob[cod_anc,cod_prob,4]==arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!="")
									num_S++;
							if(arr_prob[cod_anc,cod_prob,4]!=arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!="")
									num_N++;
							cof=(arr_prob[cod_anc,cod_prob,4]!=arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!="") ? arr_prob[cod_anc,cod_prob,3] : 1;

							sum_S1=(arr_prob[cod_anc,cod_prob,4]==arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!=""&&arr_prob[cod_anc,cod_prob,3]==0) ? arr_prob_1[cod_anc,1] : (arr_prob[cod_anc,cod_prob,4]!=arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!="") ? arr_prob[cod_anc,cod_prob,1] : 0;
							sum_N1=(arr_prob[cod_anc,cod_prob,4]==arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!=""&&arr_prob[cod_anc,cod_prob,3]==0) ? arr_prob_1[cod_anc,2] : (arr_prob[cod_anc,cod_prob,4]!=arr_prob[cod_anc,cod_prob,5]&&arr_prob[cod_anc,cod_prob,4]!="") ? arr_prob[cod_anc,cod_prob,2] : 0;
							sum_S+=sum_S1/cof;
							sum_N+=sum_N1/cof;
							sum_all+=(sum_S1+sum_N1>0) ? 1 : 0;
#							print j, cod_anc, cod_prob, arr_prob[cod_anc,cod_prob,4], arr_prob[cod_anc,cod_prob,5], sum_S1, sum_N1, cof, num_S, num_N;
						}
						print gene "\t" i "\t" arr_name[i] "\t" gene_name "\t" num_S "\t" num_N "\t" sum_S "\t" sum_N "\t" num_all "\t" sum_all;
						
					}
				}
			}
			gene=$2;
			pos=0;
		}
		pos++;
		arr_seq[pos]=$10;
		arr_name[pos]=$3;
		if($3==name)
			num_anc=pos;
		if($3==name_ref)
		{
			num_ref=pos;
			gene_name=$5;
		}
	}
}
