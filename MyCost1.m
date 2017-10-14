function z=mainFitness(x)        
    global nTask;
    global nResource;
    
     [m Secuence]=sort(x);
     SeqRes=[];
     for i=1:nTask
         SeqRes(i)=mod(Secuence(i),nResource)+1;
     end
          
      z1=Makespan(SeqRes);
      z2=Relaiability(SeqRes);
      
      
   z=[z1
      z2];  


end