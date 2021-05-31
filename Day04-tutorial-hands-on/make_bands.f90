PROGRAM  MAKE_BANDS
       character(len=1) char1, char2
       double precision  Xk_pred(3),Xk_run(3)
       double precision, dimension(:,:), allocatable  ::band
       double precision, dimension(:),     allocatable  ::Sk     
       integer nK,nband,ncol   
       character(len=80) FMT
       character(len=3) number      
          
       read(5,*) nband,nK
       ncol=nband 
                         
       allocate(band(nK,nband),Sk(nK)) 
 
        Sk(1)=0.0
        do ik=1,nK    	   
	    read(5,*) Xk_run(1),Xk_run(2),Xk_run(3)	    
            read(5,*) (band(ik,ne),ne=1,nband)	
             if(ik>1) then
	     Sk(ik)=Sk(ik-1)+sqrt( (Xk_run(1)-Xk_pred(1))**2+(Xk_run(2)-Xk_pred(2))**2+(Xk_run(3)-Xk_pred(3))**2)  
             endif
             Xk_pred(1)=Xk_run(1)
             Xk_pred(2)=Xk_run(2)
             Xk_pred(3)=Xk_run(3)    
       end do           
       
       
16     format(I3)     
 
      	write(number,16) ncol+1		
	FMT=TRIM('('// number// 'F8.3' // ')' )	  	                   
       
         write(6,*) "@# k   band (eV)"       
          do ik=1,nK           
           write(6,FMT) Sk(ik),(band(ik,ne),ne=1,ncol)                
      end do       
       
       
      
       
       
       
       deallocate(band,Sk)
       
               

         END
