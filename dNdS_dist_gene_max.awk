#################################################
# dNdS distribution
# gene maximal ratio [dN/dS]
#################################################
{
	pos++;
	arr[pos,1]=$1;
	arr[pos,2]=$5;
	arr[pos,3]=$6;
	arr[pos,4]=$7;
	arr[pos,5]=$8;
	arr[pos,6]=$10;
	arr[pos,7]=$9;
	arr[pos,8]=$3;
	arr[pos,9]=$4;
}
END{
	gene=arr[1,1];
	printf("#ID\tDs\tDn\tDs_av\tDn_av\tNsyn\tSyn\tLength\tDn/Ds\tDn/Dsav\tAccur\tgene\tStrain\n");

	for(i=1;i<=pos;i++)
	{
		if(gene!=arr[i,1])
		{
			printf("%04d\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%d\t%d\t%d\t%7.5f\t%7.5f\t%7.5f\t%s\t%s\n", gene, sum_S, sum_N, sum_S_av/num, sum_N_av/num, num_S, num_N, num_all, ratio_max, ratio_av/num, accuracy, name, strain);
			sum_S=sum_N=ratio=ratio_max=score=score_max=0;
			sum_S_av=sum_N_av=ratio_av=num=num_S=num_N=0;
			gene=arr[i,1]
		}
		ratio=(arr[i,2]!=0&&arr[i,3]!=0) ? (arr[i,3]/arr[i,5])/(arr[i,2]/arr[i,4]) : 0;
		score=arr[i,2]+arr[i,3];
		ratio_av+=ratio;
		sum_S_av+=arr[i,2]/(arr[i,4]*(arr[i,4]+arr[i,5])/arr[i,6]);
		sum_N_av+=arr[i,3]/(arr[i,5]*(arr[i,4]+arr[i,5])/arr[i,6]);
		num++;
		if(ratio>ratio_max||ratio==ratio_max&&score>score_max)
		{
			sum_S=arr[i,2]/(arr[i,4]*(arr[i,4]+arr[i,5])/arr[i,6]);
			sum_N=arr[i,3]/(arr[i,5]*(arr[i,4]+arr[i,5])/arr[i,6]);
			num_S=arr[i,2];
			num_N=arr[i,3];
			num_all=arr[i,6];
			ratio_max=ratio;
			accuracy=arr[i,6]/arr[i,7];
			strain=arr[i,8];
			name=arr[i,9];
		}
		
	}
	printf("%04d\t%7.5f\t%7.5f\t%7.5f\t%7.5f\t%d\t%d\t%d\t%7.5f\t%7.5f\t%7.5f\t%s\t%s\n", gene, sum_S, sum_N, sum_S_av/num, sum_N_av/num, num_S, num_N, num_all, ratio_max, ratio_av/num, accuracy, name, strain);
}