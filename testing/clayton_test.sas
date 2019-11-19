/*******************************************************************************
| SAS Version:  9.4
| Created By:   Thomas Drury
| Description:  Tests Clayton copula FCMP functions in copula_tools 
********************************************************************************/;

options source2;
%include "C:\Users\tad66240\OneDrive - GSK\statistics\repositories\sas\gtl_tools\contourplot1\contourplot1.sas";
%include "C:\Users\tad66240\OneDrive - GSK\statistics\repositories\sas\gtl_tools\surfaceplot1\surfaceplot1.sas";
%include "C:\Users\tad66240\OneDrive - GSK\statistics\repositories\sas\copula_tools\copula_tools.sas";
/**/
/**** INCLUDE TOOLS ***;*/
/*%include_url(url=%str(https://raw.githubusercontent.com/squiffy-statto/copula_tools/master/copula_tools.sas));*/


data clayton_test1;
  
  array u[2];
  theta = 4;

  do i = 1 to 2000;
    call clayton_bvc_sim(u1, u2, theta);
    z1 = quantile("normal", u1, 0, 1);
    z2 = quantile("normal", u2, 0, 1);
    output;
  end;

run;

proc sgplot data = clayton_test1;
  scatter x=z1 y=z2;
run;


data clayton_test2;
 
  theta = 3;
  do u1 = 0.01 to 0.99 by 0.01;
  do u2 = 0.01 to 0.99 by 0.01;
    copula_cdf = clayton_bvc_cdf(u1,u2,theta);  
    copula_pdf = clayton_bvc_pdf_cc(u1,u2,theta);  
    if u2 = 0.5 then copula_hdf = clayton_bvc_hdf_c(u1,u2,theta);  
	joint_pdf = copula_pdf*pdf("normal",probit(u1),0,1)*pdf("normal",probit(u2),0,1);
    output;
  end;
  end;

run;


proc sgrender data=clayton_test2 template=SurfaceTmplt; 
   dynamic _X='u1' _Y='u2' _Z='copula_pdf';
run;

proc sgrender data=clayton_test2 template=ContourPlotParm;
  dynamic _X="u1" _Y="u2" _Z="copula_cdf";
run;

proc sgrender data=clayton_test2 template=ContourPlotParm;
  dynamic _X="u1" _Y="u2" _Z="copula_pdf";
run;

proc sgrender data=clayton_test2 template=ContourPlotParm;
  dynamic _X="u1" _Y="u2" _Z="joint_pdf";
run;
