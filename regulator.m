function x = regulator(H,g,l,u,A,bl,bu,xinit)
     
    [x,fval,exitflag,output,lambda] = quadprog(H,g,[A; -A],[bu; -bl],[],[],l,u,xinit,options)

end




