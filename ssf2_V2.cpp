#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List ssfcpp(List x)

{

  //Parameter and variable declaration--------------------------------------------------------------------
  double imax = x["imax"]; //Iteration times
  double tmax = x["tmax"]; //Projection period
  double recfun = x["recfun"]; //Recruitment function 1:normal, 2:Autocorrelation
  double manage = x["manage"]; //Management scenario 1:Constant F, 2:Constant catch
  double catb = x["catb"]; //Average catch
  double upto = 0.0001; //Difference between estimated catch and target catch
  double thresh = 0; //Threshold of the target catch
  double age = x["age"]; //Maximum age (start age 0)
  double gender = x["gender"]; //Number of gender
  double rho = x["rho"]; //Number of gender
  NumericVector faa = x["f"]; //Fishing mortality at age
  NumericVector faa2 ((age+1)*gender); //Fishing mortality at age under constant catch management
  NumericVector maa = x["maa"]; //Natural mortarity at age
  NumericVector zaa = ((age+1)*gender); //Total mortarity at age
  NumericVector int_n = x["int_n"]; //Initial number at age
  NumericVector waa = x["waa"]; //Wright at age
  NumericVector mataa = x["mat"]; //Maturity at age
  NumericMatrix sb (imax, tmax+1); //Spawning biomass
  NumericMatrix tb (imax, tmax+1); //Total biomass
  NumericMatrix R (imax, tmax+1); //Number of recruitment fish
  NumericMatrix C (imax, tmax+1); //Catch weight
  NumericMatrix C_est (imax, tmax+1); //Estimated catch weigt under constant catch scenario
  NumericMatrix N ((age+1)*gender, tmax+1); //Population number
  NumericMatrix Sigma_R(imax, tmax+1); //Uncertainty of recruitment
  NumericMatrix n(imax, tmax+1); //Uncertainty of recruitment with enviromental effect
  NumericMatrix fmulti (imax, tmax+1); //F multiplier to operate constant catch
  double Rzero = x["Rzero"]; //R zero
  double SBzero = x["SBzero"]; //SB zero
  double h = x["h"]; //Steepness
  double sigmar = x["sigmar"]; //Sigma R
  NumericVector tmp_N ((age+1)*gender); //Temporary population number in the end of year

  //Set the initial popoulation number
  N(_,0) = clone(int_n);

  //Set the recruitment deviation
  for (int i = 0; i < imax; i++) {
    Sigma_R(i, _) = rnorm(tmax+1, 0, sigmar);
  }

  //Start calculation----------------------------------------------------------------------------------
  for (int i = 0; i < imax; i++) {
    for (int t = 0; t < tmax; t++) {

      //Calculate total mortality
      if (manage == 1) {
        zaa = exp(-(faa+maa));
      } else {
        for (int x = 0; x < 3000; x++) {
          fmulti(i,t) = x*0.001;
          faa2 = fmulti(i,t)*faa;
          C_est(i,t) =  sum(((faa2/(faa2+maa))*(1-exp(-(faa2+maa)))*N(_,t))*waa);
          thresh = catb - C_est(i,t);
          if (thresh < upto) break;
        }
        zaa = exp(-(faa2+maa));
      }

    //Population dynamics
      tmp_N = N(_,t)*zaa; //Population number in the end of year

      if (gender == 1) {
        for (int j = 0; j < age-1; j++) {
          N(j+1,t+1) = tmp_N[j];
        }
        N(age,t+1) = tmp_N[age-1]+tmp_N[age]; //Calculate plus group
        } else {
          for (int j = 0; j < age-1; j++) {
            N(j+1,t+1) = tmp_N[j];
          }
          for (int j = age+1; j < age*2; j++) {
            N(j+1,t+1) = tmp_N[j];
          }
          N(age,t+1) = tmp_N[age-1]+tmp_N[age]; //Calculate Female plus group
          N(age*2+1,t+1) = tmp_N[age*2]+tmp_N[age*2+1]; //Calculate Male plus group
        }

        //Calculate spawning and total biomass
        sb(i,0) = sum(waa*mataa*N(_,0)); //Spawning Biomass in time 0
        sb(i,t+1) = sum(waa*mataa*N(_,t+1)); //Spawing Biomass in time t+1
        tb(i,t+1) = sum(waa*N(_,t+1)); //Total Biomass in time t+1

        //Calculate recruitment
        if (recfun == 1) {
          R(i,t+1) = ((4*h*Rzero*sb(i,t+1))/(SBzero*(1-h)+sb(i,t+1)*(5*h-1)))*exp(Sigma_R(i,t+1)-(pow(sigmar,2))/2); //Recruitment in time t+1
        } else {
          n(i,0) = Sigma_R(i,0); //Coefficient representing first-order autocorrelation in time 0
          n(i,t+1) = rho*n(i,t) + pow((1-pow(rho,2)),0.5)*Sigma_R(i,t+1); //Coefficient representing first-order autocorrelation in time t+1
          R(i,t+1) = ((4*h*Rzero*sb(i,t+1))/(SBzero*(1-h)+sb(i,t+1)*(5*h-1)))*exp(n(i,t+1)-(pow(sigmar,2))/2); //Recruitment with enviroment effect in time t+1
        }

        if (gender == 1) {
          N(0,t+1) = R(i,t+1); //Number of recruitment
        } else {
          N(0,t+1) = 0.5*R(i,t+1); //Number of Female recruitment
          N(16,t+1) = 0.5*R(i,t+1); //Number of Male recruitment
        }

        //Calculate total catch amount
        if (manage == 1) {
          C(i,t) = sum(waa*((faa/(faa+maa))*(1-exp(-faa-maa))*N(_,t))); //Catch weight in time t+1 of constant F secenario
        } else {
          C(i,t) = sum(waa*((faa2/(faa2+maa))*(1-exp(-faa2-maa))*N(_,t))); //Catch weight in time t+1 of constant catch secenario
        }

      }
    }

    //Output tables---------------------------------------------------------------------------------------------------
    return List::create(_["faa"]=faa, _["maa"]=maa, _["zaa"]=zaa, _["int_n"]=int_n, _["waa"]=waa, _["mataa"]=mataa,
    _["N"]=N, _["R0"]=Rzero, _["SB0"]=SBzero, _["Steepness"]=h, _["SB"]=sb, _["TB"]=tb, _["R"]=R, _["tmp_N"]=tmp_N,
    _["C"]=C, _["Sigma_R"]=Sigma_R, _["n"]=n, _["C_est"]=C_est, _["fmulti"]=fmulti);

}
