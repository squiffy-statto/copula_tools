/*******************************************************************************
| SAS Version:  9.4
| Created By:   Thomas Drury
| Description:  Creates FCMP Functions to calculate PDFs, CDFs
|               and simulate values from different copulas.
|               These can be used for joint simulation and likelihood 
|               creation. 
|--------------------------------------------------------------------------------
| FCMP Functions List:
|--------------------------------------------------------------------------------
| Name     : gauss_bvc_cdf(u,v,o) 
| Purpose  : CDF for Gaussian bivariate copula 
| Arguments: u [REQ] = First uniform margin
|            v [REQ] = Second uniform margin 
|            o [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_bvc_hdf_c(u,v,o) 
| Purpose  : Conditional CDF given continous v for Gaussian bivariate copula 
| Arguments: u [REQ] = First uniform margin
|            v [REQ] = Second uniform margin 
|            o [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_bvc_pdf_cc(u,v,o) 
| Purpose  : PDF for continous u and v for Gaussian bivariate copula 
| Arguments: u [REQ] = First uniform margin 
|            v [REQ] = Second uniform margin
|            o [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_bvc_pdf_dd(upos,uneg,vpos,vneg,o) 
| Purpose  : PDF for discrete u and v for Gaussian bivariate copula 
| Arguments: upos [REQ] = Actual value of first margin (e.g. upos = CDF(x))
|            uneg [REQ] = Previous value of first margin (e.g. uneg = CDF(x-1))
|            vpos [REQ] = Actual value of second margin (e.g. vpos = CDF(y))
|            vneg [REQ] = Previous value of second margin (e.g. vneg = CDF(y-1))
|            o    [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_bvc_pdf_dc(upos,uneg,v,o) 
| Purpose  : PDF for discrete u and continous v for Gaussian bivariate copula 
| Arguments: upos [REQ] = Actual value of first margin (e.g. upos = CDF(x))
|            uneg [REQ] = Previous value of first margin (e.g. uneg = CDF(x-1))
|            v    [REQ] = Second uniform margin
|            o    [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_bvc_pdf_cd(u,vpos,vneg,o) 
| Purpose  : PDF for continous u and discrete v for Gaussian bivariate copula 
| Arguments: u    [REQ] = First uniform margin
|            vpos [REQ] = Actual value of first margin (e.g. upos = CDF(x))
|            vneg [REQ] = Previous value of first margin (e.g. uneg = CDF(x-1))
|            o    [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_bvc_sim(u,v,o) 
| Purpose  : Simulates pairs of u and v from Gaussian bivariate copula 
| Arguments: u [REQ] = Variable to hold first margin (must be initialized)
|            v [REQ] = Variable to hold second margin (must be initialized)
|            o [REQ] = Gauss copula dependence parameter  
|
|---------------------------------------------------------------------------------
| Name     : gauss_mvc_pdf_c(u,o) 
| Purpose  : PDF for continous vector u for Gaussian n-variate copula 
| Arguments: u [REQ] = 1 Dim array with mv copula variables 
|            o [REQ] = 2 Dim array with Gauss copula dependence matrix  
|
|---------------------------------------------------------------------------------
| Name     : gauss_mvc_sim(u,o) 
| Purpose  : Simulates multivariate values from Gaussian copula 
| Arguments: u [REQ] = 1 Dim array with mv copula variables
|            o [REQ] = 2 Dim array with Gauss copula dependence matrix  
|
***********************************************************************************/;

proc fcmp outlib = work.functions.copulas;


  **********************************************************************************;
  *** FCMP GAUSS BIVARIARIATE COPULA FUNCTIONS                                   ***;
  **********************************************************************************;
 
  *** GAUSS COPULA CDF ***;
  function gauss_bvc_cdf(u,v,o);
    if (u=0) or (v=0) then _cdf = 0;                *** COPULA PROPERTY 1:  C(0,V) = C(U,0) = C(0,0) = 0 ***;
    else if (u=1) and (v=1) then _cdf = 1;          *** COPULA PROPERTY 2:  C(1,1) = 1 ***;
    else if (0<u<1) and (v=1) then _cdf = u;        *** COPULA PROPERTY 3:  C(U,1) = U ***;
    else if (u=1) and (0<v<1) then _cdf = v;        *** COPULA PROPERTY 4:  C(1,V) = V ***;
    else if (0<u<1) and (0<v<1) then do;            *** IF BETWEEN (0,1) FOR U AND V   ***;
      _qu  = probit(min(max(u,1.0e-10),1-1.0e-10));
      _qv  = probit(min(max(v,1.0e-10),1-1.0e-10));
      _cdf = probbnrm(_qu,_qv,o);     
    end;
    return (_cdf);
  endsub;

  
  *** GAUSS COPULA PDF FOR CONTINUOUS MARGINS ***;
  function gauss_bvc_pdf_cc(u,v,o);
    _qu  = probit(min(max(u,1.0e-10),1-1.0e-10));
    _qv  = probit(min(max(v,1.0e-10),1-1.0e-10));
    _pdf = (1/sqrt(1-(o**2))) * 
           exp(-0.5 * (1/(1-o**2)) * (_qu**2 + _qv**2 - 2*o*_qu*_qv)) * 
           exp(0.5 * (_qu**2 + _qv**2));
    return(_pdf);
  endsub;  


  *** CONDITIONAL GAUSS COPULA CDF (CALLED AN H-FUNCTION) WHEN CONDITIONING VAR IS CONTINUOUS ***;
  function gauss_bvc_hdf_c(u,v,o);
    if u=0 then _hdf = 0;                       
    else if u=1 then _hdf = 1;      
    else if (0<u<1) then do;    
     _qu  = probit(min(max(u,1.0e-10),1-1.0e-10));
     _qv  = probit(min(max(v,1.0e-10),1-1.0e-10));
     _hdf = probnorm((_qu - o*_qv) / sqrt(1-o**2));
    end;
    return(_hdf);
  endsub;


  *** GAUSS COPULA PROBABILITY MASS FUNCTION FOR DISCRETE-DISCRETE VARS ***;
  function gauss_bvc_pdf_dd(upos,uneg,vpos,vneg,o);
    _pdf = ( gauss_bvc_cdf(upos,vpos,o) - 
             gauss_bvc_cdf(upos,vneg,o) - 
             gauss_bvc_cdf(uneg,vpos,o) + 
             gauss_bvc_cdf(uneg,vneg,o) ) / 
             ((upos-uneg)*(vpos-vneg)); 
    return(_pdf);
  endsub;


  *** GAUSS COPULA PROBABILITY DENSITY/MASS FUNCTION FOR DISCRETE-CONTINUOUS CASE ***;
  function gauss_bvc_pdf_dc(upos,uneg,v,o);
     _pdf = (gauss_bvc_hdf_c(upos,v,o) - gauss_bvc_hdf_c(uneg,v,o)) / (upos-uneg);
     return(_pdf);
  endsub;


  *** GAUSS COPULA PROBABILITY DENSITY/MASS FUNCTION FOR CONTINUOUS-DISCRETE VARS - SWAP ARGUMENTS IN THE GAUSS_BVC_HDF_C FUNCTION ***;
  function gauss_bvc_pdf_cd(u,vpos,vneg,o);
     _pdf = (gauss_bvc_hdf_c(vpos,u,o) - gauss_bvc_hdf_c(vneg,u,o)) / (vpos-vneg);
     return(_pdf);
  endsub;


  *** SIMULATE BIVARIATE GAUSS COPULA PAIR ***;
  subroutine gauss_bvc_sim(u,v,o);
    outargs u,v;
    _z1 = rand("normal",0,1);
    _z2 = rand("normal",0,1);
    u = probnorm(_z1);
    v = probnorm(sqrt(1-o**2)*_z2 + (o*_z1));
  endsub;

  **********************************************************************************;
  *** FCMP GAUSS N-VARIATE COPULA FUNCTIONS                                      ***;
  **********************************************************************************;
  
  *** N-VARIATE GAUSS COPULA PDF FOR CONTINUOUS MARGINS  ***;
  function gauss_mvc_pdf_c(u[*],o[*,*]);
   udim1 = dim1(u); 
   odim1 = dim1(o);
   odim2 = dim2(o);
   if not ( udim1 = odim1 = odim2 ) then do;
     msg1 = "ER"||upcase("ror:(FCMP):")||"The Function GAUSS_MVC_PDF_C does not have matching array sizes.";
     msg2 = "ER"||upcase("ror:(FCMP): Dimensions:");
     put msg1;
     put msg2 udim1= odim1= odim2=;
     _val = .e;
   end;
   else do;
     n = udim1;
     array qvec[1,1] / nosymbols;     
     array omat[1,1] / nosymbols;
     array qtrs[1,1] / nosymbols;
     array oinv[1,1] / nosymbols;
     array iden[1,1] / nosymbols;
     array mat1[1,1] / nosymbols;
     array vec1[1,1] / nosymbols;
     array val1[1,1] / nosymbols;
     call dynamic_array(qvec,n,1);
     call dynamic_array(omat,n,n);
     call dynamic_array(qtrs,1,n);
     call dynamic_array(oinv,n,n);     
     call dynamic_array(iden,n,n);
     call dynamic_array(mat1,n,n);     
     call dynamic_array(vec1,n,1);
     do ii = 1 to n;
       qvec[ii,1] = probit(min(max(u[ii],10e-10),1-10e-10));
       qtrs[1,ii] = qvec[ii,1];
       do jj = 1 to n;
         omat[ii,jj] = o[ii,jj];
       end;
     end;
     call det(omat,deto);
     call identity(iden);
     call inv(omat,oinv);
     call subtractmatrix(oinv,iden,mat1);
     call mult(mat1,qvec,vec1);
     call mult(qtrs,vec1,val1);
     _val = (1/sqrt(deto))*exp(-0.5*val1[1,1]);
   end;
   return(_val);
  endsub; 


  *** SIMULATE VALUES FROM N-DIM GAUSS COPULA ***;
  subroutine gauss_mvc_sim(u[*],o[*,*]);   
    outargs u;
    udim1 = dim1(u); 
    odim1 = dim1(o);
    odim2 = dim2(o);
    if not ( udim1 = odim1 = odim2 ) then do;
      msg1 = "ER"||upcase("ror:(FCMP):")||"The Function GAUSS_MVC_SIM does not have matching array sizes.";
      msg2 = "ER"||upcase("ror:(FCMP): Dimensions:");
      put msg1;
      put msg2 udim1= odim1= odim2=;
    end;
    else do;
     n = udim1;
     array uvec[1,1] / nosymbols;     
     array omat[1,1] / nosymbols;
     array chol[1,1] / nosymbols;
     array cvec[1,1] / nosymbols;
     call dynamic_array(uvec,n,1);     
     call dynamic_array(omat,n,n);
     call dynamic_array(chol,n,n);     
     call dynamic_array(cvec,n,1);
     do ii = 1 to n;
       uvec[ii,1] = rand("NORMAL",0,1);
       do jj = 1 to n;
         omat[ii,jj] = o[ii,jj];
       end;
     end;
     call chol(omat,chol);
     call mult(chol,uvec,cvec);
     do ii = 1 to n;
       u[ii] = cdf("normal", cvec[ii,1], 0, 1); 
     end;
    end;
  endsub;


  **********************************************************************************;
  *** FCMP BIVARIATE T COPULA FUNCTIONS                                          ***;
  **********************************************************************************;
   
  *** TO ADD ***;

run;
quit;

options cmplib = work.functions;
