library(pryr)
n<-200
rp<-1000 #循环次数
n_var<-10
set.seed(123)
r<-runif(n,min=0,max=1)
f_creatdata<-function(n_var){
  ind<-rep(1:n)
  d_pat<-data.frame(ind)
  
  for(i in 1:n_var){
    x<-runif(n,min=0,max=1)
    y<-(x>0.5)
    z<-y+1
    d_pat<-cbind(d_pat,z)
    names(d_pat)[i+1]<-paste("X",i,sep="") 
  }
  
  d_pat$Group1<-""
  d_pat$Group2<-""  
  d_pat$Group3<-""
  return(d_pat)
}
control_limit <- 0.3
vote_f <- function(table, expected, actual) {
  r <- runif(1)
  p_value <- chisq.test(table)$p.value
  if (p_value > control_limit) {
    if (r < 0.5) {
      return("Vote for A")
    } else {
      return("Vote for B")
    }
  } else {
    if (expected[1] > actual[1]) {
      return("Vote for A")
    } else if (expected[2] > actual[2]) {
      return("Vote for B")
    } else {
      if (r< 0.5) {
        return("Vote for A")
      } else {
        return("Vote for B")
      }
    }
  }
}
vote_result <- function(votes, bias_prob = 0.65) {
  votes_for_A <- sum(votes == "Vote for A")
  votes_against_A <- sum(votes == "Vote for B")
  if (votes_for_A > votes_against_A) {
    return(ifelse(runif(1) < bias_prob, "A", "B"))
  } else if (votes_against_A > votes_for_A) {
    return(ifelse(runif(1) < bias_prob, "B", "A"))
  } else {
    return(sample(c("A", "B"), 1, replace = TRUE))
  }
}
new_patient_allocation <- function(d_pat,n_var) {
  votes <- c()
  for(i in 1:n_var){
    var_name<-paste("X",i,sep = "")  
    table <- table(d_pat[[var_name]], d_pat$Group1)
    expected <- chisq.test(table)$expected
    actual <- table
    votes <- c(votes, vote_f(table, expected, actual))
    p_value <- chisq.test(table)$p.value  
    p_values <- c(p_values, p_value)      
  }
  return(vote_result(votes))
}
f_distribute<-function(ind,d_pat,n_var){
  dA<-d_pat[d_pat$ind<=ind,]
  dA[dA$ind==ind,"Group2"]<-"A"
  dB<-d_pat[d_pat$ind<=ind,]
  dB[dB$ind==ind,"Group2"]<-"B"
  
  cha_A<-0
  cha_B<-0
  
  for(i in 1:n_var){
    var_name<-paste("X",i,sep = "")
    cha_A<-cha_A+f_sub_chisq(var_name,2,dA)
    cha_B<-cha_B+f_sub_chisq(var_name,2,dB)
  }
  
  rm(dA)
  rm(dB)
  rm(d_pat)
  
  res<-c(cha_A,cha_B)
  return(res)
}

f_sub_chisq<-function(var_name,var_count,d){
  cell_zero<-data.frame(Var1=rep(1:var_count,each=2),Var2=rep(c("A","B")),Freq=rep(0))
  cell_count<-data.frame(table(d[[var_name]],d$Group2))
  cell_merge<-merge(cell_zero,cell_count,by=c("Var1","Var2"),all=T)
  cell_merge[is.na(cell_merge$Freq.y),"Freq.y"]<-0
  cell_merge$Freq<-cell_merge$Freq.x+cell_merge$Freq.y

  cell_expect<-sum(cell_merge$Freq)/6
  cha_phase<-sum((cell_merge$Freq-cell_expect)^2/cell_expect)

  return(cha_phase)
}
f_unbalance<-function(ind,d_pat,n_var){
  dA<-d_pat[d_pat$ind<=ind,]
  dA[dA$ind==ind,"Group3"]<-"A"
  dB<-d_pat[d_pat$ind<=ind,]
  dB[dB$ind==ind,"Group3"]<-"B"
  
  diff_A<-0
  diff_B<-0
  
  for(i in 1:n_var){
    var_name<-paste("X",i,sep = "")
    var_value<-d_pat[d_pat$ind==ind,var_name]
    diff_A<-diff_A+abs(nrow(dA[dA$Group3=="A" & dA[,var_name]==var_value,])
                       -nrow(dA[dA$Group3=="B" & dA[,var_name]==var_value,]))
    diff_B<-diff_B+abs(nrow(dB[dB$Group3=="A" & dB[,var_name]==var_value,])
                       -nrow(dB[dB$Group3=="B" & dB[,var_name]==var_value,]))
    #diff_A<-diff_A+abs(nrow(dA[dA$Group1=="A" ,])-nrow(dA[dA$Group1=="B" ,]))
    #diff_B<-diff_B+abs(nrow(dB[dB$Group1=="A" ,])-nrow(dB[dB$Group1=="B" ,]))
    
  }
      res<-c(diff_A,diff_B)
  return(res)
}

f_main<-function(){
  
  for(n_var in 1:10){
    print(paste("variable number =",n_var,":"))
    pv1<-as.data.frame(matrix(,rp,n_var+1))
    pv2<-as.data.frame(matrix(,rp,n_var+1))
    pv3<-as.data.frame(matrix(,rp,n_var+1))
    colnames(pv1)[1]<-"count_diff"
    colnames(pv2)[1]<-"count_diff"
    colnames(pv3)[1]<-"count_diff"
    
    for(i in 1:n_var){
      colnames(pv1)[i+1]<-paste("X",i,sep = "")
      colnames(pv2)[i+1]<-paste("X",i,sep = "")
      colnames(pv3)[i+1]<-paste("X",i,sep = "")
    }
    
    #csv<-as.data.frame(matrix(,0,c(1:(n_var*6+7)),ncol=n_var*6+7)) 
    k<-0
    csv<-as.data.frame(matrix(,0,n_var*6+7))    
    colnames(csv)[1]<-"sample_size"
    colnames(csv)[2]<-"diff_00"
    colnames(csv)[3]<-"diff_05"
    colnames(csv)[4]<-"diff_50"
    colnames(csv)[5]<-"diff_95"
    colnames(csv)[6]<-"diff_100"
    for(i in 1:n_var){
      k<-(i)*5
      colnames(csv)[k+2]<-paste("X",i,"_00",sep = "")
      colnames(csv)[k+3]<-paste("X",i,"_05",sep = "")
      colnames(csv)[k+4]<-paste("X",i,"_50",sep = "")
      colnames(csv)[k+5]<-paste("X",i,"_95",sep = "")
      colnames(csv)[k+6]<-paste("X",i,"_100",sep = "")
    }
    i<-i+1
    k<-(i)*5+2
    colnames(csv)[k]<-"diff_mean"
    for(i in 1:n_var){
      colnames(csv)[k+i]<-paste("X",i,"_mean",sep = "")
    }
    
    
    
    for(j in 1:rp){
      d_pat<-f_creatdata(n_var)
      d_pat$Group1[1:20]<-sample(c("A", "B"),20, replace = TRUE)
      
      for(i in 20:(n-1)){
        d_part<-d_pat[1:i,]
        t_news<-new_patient_allocation(d_part,n_var)
        d_pat$Group1[i+1]<-t_news
      }
      
     
      for(i in 1:n){
        print(c(i,"-----------------------------"))
        if(i==1) { 
          if(r[i]<0.5) {d_pat[i,"Group2"]<-"A";d_pat[i,"Group3"]<-"A";}
          else {d_pat[i,"Group2"]<-"B";d_pat[i,"Group3"]<-"B";}
        }
        else { 
      u2<-f_distribute(i,d_pat,n_var)
          print(cat(" u2:",u2))
          
          if(u2[1]==u2[2]){
            if(r[i]<0.5) d_pat[i,"Group2"]<-"A"
            else d_pat[i,"Group2"]<-"B"
          }else if(u2[1]>u2[2]){
            if(r[i]<0.8) d_pat[i,"Group2"]<-"B"
            else d_pat[i,"Group2"]<-"A"
          }else if(u2[1]<u2[2]){
            if(r[i]<0.8) d_pat[i,"Group2"]<-"A"
            else d_pat[i,"Group2"]<-"B"
          }
          u3<-f_unbalance(i,d_pat,n_var)
          print(cat(" u3:",u3))
          
          if(u3[1]==u3[2]){
            if(r[i]<0.5) d_pat[i,"Group3"]<-"A"
            else d_pat[i,"Group3"]<-"B"
          }else if(u3[1]>u3[2]){
            if(r[i]<0.8) d_pat[i,"Group3"]<-"B"
            else d_pat[i,"Group3"]<-"A"
          }else if(u3[1]<u3[2]){
            if(r[i]<0.8) d_pat[i,"Group3"]<-"A"
            else d_pat[i,"Group3"]<-"B"
          }
        }
      }
      
      
      tmp1<-table(d_pat[,"Group1"])
      pv1[j,1]<-abs(tmp1[1]-tmp1[2])
      
      tmp2<-table(d_pat[,"Group2"])
      pv2[j,1]<-abs(tmp2[1]-tmp2[2])
      
      tmp3<-table(d_pat[,"Group3"])
      pv3[j,1]<-abs(tmp3[1]-tmp3[2]) 
      
      for(i in 1:n_var){
        var_name<-paste("X",i,sep = "")
        
        t1<-chisq.test(table(d_pat[,c(var_name,"Group1")]),correct=TRUE)
        pv1[j,i+1]<-t1$p.value 
        
        t2<-chisq.test(table(d_pat[,c(var_name,"Group2")]),correct=TRUE)
        pv2[j,i+1]<-t2$p.value
        
        t3<-chisq.test(table(d_pat[,c(var_name,"Group3")]),correct=TRUE)
        pv3[j,i+1]<-t3$p.value
      }
      
      rm(d_pat)
      
      print(paste("the",j,"th loop done."))
    }
    q1<-apply(pv1, 2, quantile,probs = c(0,0.05,0.5,0.95,1))
    m1<-apply(pv1, 2, mean)
    
    
    q2<-apply(pv2, 2, quantile,probs = c(0,0.05,0.5,0.95,1))
    m2<-apply(pv2, 2, mean)
    
    q3<-apply(pv3, 2, quantile,probs = c(0,0.05,0.5,0.95,1))
    m3<-apply(pv3, 2, mean)
    res<-rbind(c(n,matrix(q1,1,),m1),c(n,matrix(q2,1,),m2),c(n,matrix(q3,1,),m3))
    
    csv<-rbind(csv,res)
    
    write.csv(pv1,file=paste("pv1_",n_var,".csv"))
    write.csv(pv2,file=paste("pv2_",n_var,".csv"))
    write.csv(pv3,file=paste("pv3_",n_var,".csv"))
    write.csv(csv,file=paste("pv4_",n_var,".csv"))
    
   }
  
  
  min_m2<-min(m2)
  return(min_m2)
}

f_main()