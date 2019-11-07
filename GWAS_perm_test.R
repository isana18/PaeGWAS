## This code was used to perform GWAS tests and corresponding permutation tests to identify genes that associate with high or low virulence 


## Function takes @row in the matrix, @virulence vector and performs statistical tests for this row
testAssociationPerGene <- function(row,virulence,indexes_all)
{
    indexes_with_gene = which(row > 0)
    indexes_without_gene = indexes_all[-indexes_with_gene]
	
    with_gene_virulence = virulence[indexes_with_gene]
    without_gene_virulence = virulence[indexes_without_gene]

    if(length(with_gene_virulence)>0 && length(without_gene_virulence)>0) 
    { 
        ## wilcoxon test
    	pval1_1 = wilcox.test(with_gene_virulence,without_gene_virulence,alternative = "less")$p.value
	pval1_2 = wilcox.test(with_gene_virulence,without_gene_virulence,alternative = "greater")$p.value
    
        scores = virulence
        labels = row
        labels[which(labels > 0)] = 1
        ## linear regression
        fit = lm(scores~labels)
        pval_2 = summary(fit)$coefficients[2,4]
        slope = summary(fit)$coefficients[2,1] 

    }else {
	pval1_1 = 1
	pval1_2 = 1
        pval_2=1
	slope = 0
    }

    if(slope<0)
    {	
       log_LR = -log10(pval_2)
    }else {
       log_LR = log10(pval_2)
    }

    if(pval1_1<pval1_2)
    {
       log_MW = -log10(pval1_1)
    }else {
       log_MW = log10(pval1_2)
    }

    return(c(log_MW,log_LR))
}

## Takes @matrix of clusters*strains, @virulence vector and @indexes of the strains, and computes statistical associations for the whole matrix 
testAssociation <- function(matrix,virulence,indexes_all)
{
	res = c()
	for(i in 1:dim(matrix)[1])
	{
		vec = testAssociationPerGene(matrix[i,],virulence,indexes_all)
		if(i==1)
		{
		   res = vec
		}
		else
		{
		   res = rbind(res,vec)	
		}
	}
        colnames(res) = c("logMW","logLR")
	return (res)
}



####################  PERMUTATION TEST on the gene level ##############################33

updateCounts <- function(originalPval,shuffPval,counts)
{
	length = dim(originalPval)[1]
	for (i in 1:length){
		#print(i)
		if((originalPval[i,1]<0 && originalPval[i,1]>shuffPval[i,1]) || (originalPval[i,1]>0 && originalPval[i,1]<shuffPval[i,1])){
                    counts[i,1] = counts[i,1]+1
		    #print("true")	
		}
	        if((originalPval[i,2]<0 && originalPval[i,2]>shuffPval[i,2]) || (originalPval[i,2]>0 && originalPval[i,2]<shuffPval[i,2])){
                    counts[i,2] = counts[i,2]+1
		    #print("true")
		}
	}
        return (counts)
}

performPermutationTest <-function(matrix,virulence,indexes_all, originalPvals, R)
{ 
	counts = matrix(0,dim(originalPvals)[1],2)
	for (i in 1:R) {
		print(i)
		virulence_sampled = sample(virulence, size=length(virulence), replace=F)
		res_shuff = testAssociation(matrix,virulence_sampled,indexes_all)
		counts = updateCounts(originalPvals,res_shuff,counts)
	}
        
	pvals = counts/R 
  	colnames(pvals) = c("pval_logMW","pval_logLR")
        return(pvals)
}



## Read presence/absence file
inputFile = "InputMatrix.txt"
data =  read.csv(inputFile,sep ="\t", header = TRUE)
data_x = data.matrix(data[,2:dim(data)[2]])
rownames(data_x) = data[,1]

## virulence info
virulence = c(3.354,2.465,4.817,2.730,1.462,5.035,8.725,3.900,4.628,8.210,6.328,3.235,3.472,5.378,8.270,2.691,2.522,10.495,6.117,3.521,6.461,9.437,5.782,6.651,5.080,9.110,6.419,2.285,4.397,5.261,3.195,4.113,5.591,9.393,10.050,5.540,3.822,4.524,7.241,1.691,2.839,5.847,2.730,5.310,3.134,6.103,8.513,2.669,6.079,8.786,2.255,6.239)

indexes_all = seq(1,52,by=1)


result = testAssociation(data_x,virulence,indexes_all)
rownames(result) = rownames(data_x)
write.csv(result,"StatisticalTests_pvals.csv")

p_vals = performPermutationTest(data_x,virulence,indexes_all,result, 10)
results_withPermTest = cbind(result,p_vals)	
write.csv(results_withPermTest,"StatTest_withPerm_byGene_10.csv")
	
