function PlotTet(T,X)
figure 
tetramesh(T,X,'FaceAlpha',0.4)
hold on 

   N=unique(T); % all nodes
   for i=1:length(N)
       text(X(N(i),1),X(N(i),2) ,X(N(i),3),num2str(N(i)),'Color','red')
   end 

end 