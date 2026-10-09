"""
Joseph D. Taylor
Department of Physics 
University of Bath, United Kingdom
j.taylor@bath.ac.uk

June 7th 2018

Generation of the C++ files required by IPOPT.

"""

import os
import re
import sympy as sym

class BuildIpoptCpp(object):
    
    def __init__(self, DS, name='Test', specs_path='problem_info.txt'):
        """
        """
        if not os.path.exists("../outputMSS"):
            os.makedirs("../outputMSS")
        
        self.prob = name
        self.prob_path = specs_path
        self.n_vars = len(DS.sVars) + len(DS.sCtrl)
        self.n_pars = len(DS.sPars)
        self.n_data = len(DS.obs)
        self.n_constr = len(DS.dEqns)
        self.n_ctrl = len(DS.idx_ctrl)
        self.n_svar = len(DS.sVars)
        self.con_limit = 2
        self.obj_cpp = self.convert_objective(DS, DS.obj_func)
        self.build_hpp(DS)
        self.build_main(DS)
        self.build_nlp(DS)
        self.build_make( "/opt/IPOPT/lib" )
        
        
    def convert_objective(self, DS, string):
        dict_obj = {}
        joint_str = DS.sVars+DS.sCtrl
        for var in range(len(joint_str)):
            if var > 0:
                dict_obj.update({joint_str[var] : 'x[i+%d*nObs_]' % var})
            else:
                dict_obj.update({joint_str[var] : 'x[i]'})
        for data in range(len(DS.obs)):
            if DS.idx_meas[data] > 0:
                dict_obj.update({'%s' % DS.obs[data][0] : 'data%d_[i]' % (DS.idx_meas[data]+1)})
            else: 
                dict_obj.update({'%s' % DS.obs[data][0] : 'data_[i]'})
                
        dict_obj.update({'h_' : 'h[j]'})
                
        output = re.sub(r'\b' + '|'.join(dict_obj.keys()) + r'\b', 
                        lambda m: dict_obj[m.group(0)], string)
        
        return output


    def convert_constraint(self, DS, string):
        dict_con = {}
        joint_str = DS.sVars+DS.sCtrl
        for var in range(len(joint_str)):
            if var > 0:
                dict_con.update({joint_str[var] + '_t0' : 'x[j+%d*nObs_]' % var})
                dict_con.update({joint_str[var] + '_t1' : 'x[j+1+%d*nObs_]' % var})
                dict_con.update({joint_str[var] + '_t2' : 'x[j+2+%d*nObs_]' % var})
                dict_con.update({joint_str[var] + '_t3' : 'x[j+3+%d*nObs_]' % var})
                dict_con.update({joint_str[var] + '_t4' : 'x[j+4+%d*nObs_]' % var})
            else:
                dict_con.update({joint_str[var] + '_t0' : 'x[j]'})
                dict_con.update({joint_str[var] + '_t1' : 'x[j+1]'})
                dict_con.update({joint_str[var] + '_t2' : 'x[j+2]'})
                dict_con.update({joint_str[var] + '_t3' : 'x[j+3]'})
                dict_con.update({joint_str[var] + '_t4' : 'x[j+4]'})
        for par in range(len(DS.sPars)):
            if par > 0:
                dict_con.update({DS.sPars[par] : 'x[%d*nObs_+%d]' % (self.n_vars,par)})
            else:
                dict_con.update({DS.sPars[par] : 'x[%d*nObs_]' % (self.n_vars)})
        for data in range(len(DS.obs)):
            if DS.idx_meas[data] > 0:
                dict_con.update({'%s' % DS.obs[data][0] : 'data%d_[j]' % (DS.idx_meas[data]+1)})
            else:
                dict_con.update({'%s' % DS.obs[data][0] : 'data_[j]'})
        try:
            for cur in range(len(DS.forc)):
                if cur > 0:
                    dict_con.update({'%s' % (DS.forc[cur][0]): 'inj%d_[j]' % (cur+1)})
        except:
            pass
        
        dict_con.update({'h_' : 'h[j]'})
    
        output = re.sub(r'\b' + '|'.join(dict_con.keys()) + r'\b', 
                        lambda m: dict_con[m.group(0)], string)
        
        try:
            dict_con.update({'%s' % (DS.forc[0][0]) : 'inj_[j]'})
        except:
            pass
        
        output = re.sub(r'\b' + '|'.join(dict_con.keys()) + r'\b', 
                        lambda m: dict_con[m.group(0)], string)
        return output
    
        
    def build_hpp(self, DS):
        """
        """
        hppstring = self.prob + '_nlp.hpp'
        f = open(hppstring, 'w')
        
        f.write("""// %s_nlp.hpp""" % self.prob + """
// IPOPT code for variational inference
// Author:  Joseph Taylor
// Author 2021, 2022: Stephen A Wells
        
        
#ifndef __%s_NLP_HPP__""" % self.prob.upper() + """
#define __%s_NLP_HPP__""" % self.prob.upper() + """

#include "IpTNLP.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace Ipopt;

// This inherits from Ipopt's TNLP
class %s_NLP : public TNLP""" % self.prob + """
{
public:

  //SAW new properties added to control MSS behaviour
  int nObsVisible; //make the private nObs_ available
  int n_mss; // MSS covers n data points at a time
  double mss_multiplier; //increase n_mss by this factor on restart
  bool setToData; //set x[j] to data[j] on initialisation
  bool useInitdata; //take an initdata file?
  string initdataFile; //filename for initdata information
  double mssTrack; //track best value of objective/log(n_mss) in MSS loop

  bool* correctJacobian; //notes when jacobian entries must be zeroed for mss
  
  bool doYank; // if true (default), parameters at edge are yanked to middle
  
  int objFactor; //normalises objective function relative to constraints
  //should be between 1 and nObs_ ; nObs_ is previous behaviour

  /** constructor taking in problem data */
  %s_NLP();""" % self.prob + """

  /** default destructor */
  virtual ~%s_NLP();""" % self.prob + """

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing.
   *  
   */
   
  int nObs_;
  int nSkip;
  int adapt;
  double* h;
  string* prob; //array to read from problem_info.txt
  string* mssConf; //SAW added: array to read from mss_info.txt
  """) #end of first write statement, time to work on data arrays
        #SAW the data arrays should be written as before
        for var in range(len(DS.hess_obj)):
            if DS.hess_obj[var][0] == 0:
                f.write("""double* data_;
  """)
            else: 
                f.write("""double* data%s_;
  """ % (DS.hess_obj[var][0]+1))
        #now the inj arrays for injected data:
        try:        
            for inj in range(len(DS.forc)):
                if inj == 0:
                    f.write("""double* inj_;
  """)
                else: 
                    f.write("""double* inj%d_;
  """ % (inj+1))
        except:
            pass
        #SAW note: replace "init[n]_" arrays with saved_state
        f.write("""  //saved_state is used for init and also for restarts
  double* saved_state;
  """)
        # now the "bounds" etc. statements 
        f.write("""double* bounds;
  """)
                       
        f.write("""int n_vars;
  int n_ctrl;
  int n_pars;
  """)
                
        f.write("""   
  //@{
  //%s_NLP();
  %s_NLP(const %s_NLP&);
  %s_NLP& operator=(const %s_NLP&);""" % (self.prob, self.prob, 
  self.prob, self.prob, self.prob) + """
  //@}

  /** @name NLP data */

};


#endif
        """)
        
        f.close() #end of writing the .hpp file
        
    
    def build_main(self, DS):
        """
        """
        mainstring = self.prob + '_main.cpp'
        f = open(mainstring, 'w')

        f.write("""// %s_main.cpp""" % self.prob + """
// IPOPT code for variational inference
// Author:  Joseph Taylor
        
#include <cstdio>
#include <stdio.h>    
#include <stdlib.h>
#include <time.h>  
        
#include "IpIpoptApplication.hpp"
#include "%s_nlp.hpp" """ % self.prob + """
        
using namespace Ipopt;
        
int main(int argv, char* argc[])
{
          
  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<%s""" % self.prob + """_NLP> mynlp = new %s""" % self.prob + """_NLP(); 
        
  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  SmartPtr<IpoptApplication> app = new IpoptApplication();
        
  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-10);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "../outputMSS/ipopt.out");
        
  // Intialize the IpoptApplication and process the options
  app->Initialize();

  //int savedN = mynlp->n_mss;
  //cout << "Initial n_mss " << mynlp->n_mss << endl;      
  // Ask Ipopt to solve the problem
  ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);
        
  if (status == Solve_Succeeded) {
    printf("\\n\\n*** The problem solved!\\n");
  }
  else {
    printf("\\n\\n*** The problem FAILED!\\n");
  }

  //stop if mss_multiplier was set to 1 or less
  if ( !(mynlp->mss_multiplier > 1.0) ){
    printf("\\n Halting (no multiplier active).\\n");
    return (int) status; //stop here!
  }

  //loop, increasing n_mss every time (done by finalize)
  int nLoops = 0;
  while ( mynlp->n_mss < mynlp->nObsVisible ){
	  nLoops += 1;
	  cout << "Attempting loop " << nLoops << endl;
	  cout << "Current n_mss " << mynlp->n_mss << endl;
	  status = app->OptimizeTNLP(mynlp);
			
	    if (status == Solve_Succeeded) {
                printf("\\n\\n*** The problem solved!\\n");
              }
              else {
                printf("\\n\\n*** The problem FAILED!\\n");
              }
  }
  
  cout  << "Closing run (MSS becomes LSS)" << endl;
  cout << "Current n_mss " << mynlp->n_mss << endl;

  status = app->OptimizeTNLP(mynlp);
			
  if (status == Solve_Succeeded) {
    printf("\\n\\n*** The problem solved!\\n");
  }
  else {
    printf("\\n\\n*** The problem FAILED!\\n");
  }

  //stop after loop
  return (int) status;
}
""")       
        f.close()
    
    
    def build_nlp(self, DS):
        """
        """
        nlpstring = self.prob + '_nlp.cpp'
        f = open(nlpstring, 'w')

        f.write("""// Id: %s_nlp.cpp""" % self.prob + """
// IPOPT code for variational inference
// Author:  Joseph Taylor
        
#include "%s_nlp.hpp" """ % self.prob + """
        
#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
        
using namespace Ipopt;

// constructor
%s_NLP::%s_NLP()""" % (self.prob, self.prob) + """
{
  
  """)
        try:
            for forc in range(len(DS.forc)):
                if forc == 0:
                    f.write("""vector<double> injData;
  """)
                else: 
                    f.write("""vector<double> injData%d;
  """ % (forc+1))
        except:
            pass
        
        for var in range(len(DS.hess_obj)):
            if DS.hess_obj[var][0] < len(DS.sVars):
                if DS.hess_obj[var][0] == 0:
                    f.write("""vector<double> inputData;
  """)
                else: 
                    f.write("""vector<double> inputData%d;
  """ % ((DS.hess_obj[var][0]+1)))     

        f.write("""
  n_vars = %d;
  n_ctrl = %d;
  n_pars = %d;
  
  """ % (self.n_svar,self.n_ctrl,self.n_pars))
        
        f.write("""// read in problem information
  fstream fp;
  int count = 0;
  std::string line;
  prob = new string[100];
  
  ifstream myfile ("%s");""" % self.prob_path + """
  if (myfile.is_open()){
    while (getline(myfile,line)){
      if (line[0] != '#'){  
            prob[count] = line;
            count++;
      }
    }
    myfile.close();
  } 
  else {
    cout << "Unable to open problem info file." << endl;
    return; //SAW added failsafe kill
  }

  nObs_ = atoi(prob[0].c_str());
  nSkip = atoi(prob[1].c_str());
  adapt = atoi(prob[2].c_str());

  //SAW making the value of nObs_ available:
  nObsVisible = nObs_;

  mssTrack = 1.0e6; //fictitious very large first value for tracking best mss score
  
  """)
        #SAW introduce mss_info.txt readin here, based on problem_info logic
        #should be a block without any variables inserted...
        f.write("""//read in mss configuration information
  count = 0;
  mssConf = new std::string[100]; //just copying from prob :)

  ifstream myMSSfile ("mss_info.txt");
  if (myMSSfile.is_open()){
    while (getline(myMSSfile,line)){
      if (line[0] != '#'){  
            mssConf[count] = line;
            count++;
      }
    }
    myMSSfile.close();
  } 
  else {
    cerr << "Unable to open MSS info file." << endl;
    return; //failsafe!
  }
  
  n_mss = atoi( mssConf[0].c_str());
  mss_multiplier = atof( mssConf[1].c_str());
  setToData =  ( atoi( mssConf[2].c_str()) > 0 ); //true if 1 written here...

  int useInitFlag = atoi( mssConf[3].c_str()); //pick up 0,1,2
  
  useInitdata =  ( useInitFlag > 0 ); //true if 1+ written here...
  initdataFile = mssConf[4]; //grab from file, should be ../data/initdata.dat usually...
  doYank = true; //default behaviour, rescue yank

  objFactor = atoi( mssConf[5].c_str() );

  if (objFactor == 0 ) objFactor = nObs_;
  cout << "Objective normalise objFactor: " << objFactor << endl; //debug check
  

  //it's important that n_mss be even:
  if ( n_mss%2 > 0 ){
    n_mss += 1; //increment odd to even.
  }
  cout << "Initial n_mss: " << n_mss << endl;
  
""")
            
        f.write("""// read in state variable and parameter bounds
  count = 0;
  std::string li;
  bounds = new double[300];
  ifstream boundfile ("%s");""" % DS.bounds_path + """
  if (boundfile.is_open()){
    while (getline(boundfile,li)){
      if (li[0] != '#'){  
	     std::string input = li;
		  std::istringstream ss(input);
		  std::string token;
		  while(std::getline(ss, token, ',')) {
		    bounds[count] = atof((token).c_str());
		    count++;		
         }					
      }
    }
    boundfile.close();
  } 
  else {
    cout << "Unable to open bounds file." << endl;
    return;
  }
 
  // read in step size from problem info, or text file if adaptive
  h = new double[nObs_];        
  string filename;
  filename = prob[4];
  
  FILE *rFile0;
  rFile0 = fopen(filename.c_str(), "r");
            
  if (adapt == 1){                
    for (Index j = 0; j<nObs_; j++){
      if(fscanf(rFile0, "%lf", &h[j])!=EOF){};
    }
    fclose(rFile0); 
  }else{
    for (Index j = 0; j<nObs_; j++){
      h[j] = atof(prob[3].c_str());
    }
  }
   
""") #we've just finished the timestep part
        #next is the initialisation
        #MSS: initialisation now uses saved_state
        #and useInitdata option
        #first, the setup for saved_state using bounds
        f.write("""
  //since we have bounds etc:
  //set up the saved_state variable
  //this is initialisation if no initdata file is given
  
  saved_state = new Number[ ((n_vars+2*n_ctrl)*nObs_)+n_pars ];
  //exactly the same size as the x vector
  
  //set bounds value for state variables
  for (Index var = 0; var < n_vars; var++){
	  for (Index i=var*nObs_; i<(var+1)*nObs_; i++) {
		  saved_state[i] = bounds[3*var+2];
	  }
  }
  //set bounds value for control variables:
  for (Index ctrl = 0; ctrl < n_ctrl; ctrl+=2){
	  for (Index i=(n_vars+2*ctrl)*nObs_; i<(n_vars+2*ctrl+1)*nObs_; i++) {
		  saved_state[i] = bounds[3*(n_vars+2*ctrl)+2];
		  saved_state[i+nObs_] = bounds[3*(n_vars+2*ctrl+1)+2];
	  }
  }
  //set bounds values for parameters
  for (Index par = 0; par < n_pars; par++){	
    saved_state[(n_vars+2*n_ctrl)*nObs_+par] = bounds[3*(n_vars+2*n_ctrl+par)+2];
  }  
""")
        #next, if we are using an initdata file, read and use it:
        #191021 adding a restart-file function to this section
        f.write("""
  if ( useInitdata ){
	  cout << "Using initial data information from ";
	  cout << initdataFile << endl;

	  // read in initial state history guess
	  const int cols = n_vars;

	  int nToSkip = 0;
	  if ( useInitFlag == 1 ) nToSkip = nSkip;
	  //a flag of 2 will not skip any lines, for restart files
				
	  ifstream infile;
	  infile.open( initdataFile.c_str() );

	  if (!infile) {
		cerr << "Unable to open initial state file.\\n" << endl;
		return;   // call system to stop
	  }

	  double initFile[cols*(nObs_+nSkip)];

	  for (int i = 0; i < cols*(nObs_+nToSkip); i++) {
		if (!(infile >> initFile[i])) {
		  cerr << "Unexpected end of initial state file.\\n" << endl;
		  return;   // call system to stop
		}    
	  }
	  infile.close();
	  
	  //need correspondence of saved_state/x with init and initFile...
	  for ( Index i = 0 ; i < nObs_ ; i++ ){
		  for ( Index j = 0 ; j < n_vars ; j++ ){
			  saved_state[ i + j*nObs_ ] = initFile[ cols*(i+nToSkip) + j ];
		  }
	  }
	  
	  //end of the initdata section
  }
  int nRead = 0; //used for line counting and skipping, next
""")

        #now begins the "read in observations" part (not affected by MSS)
        for var in range(len(DS.hess_obj)):
            if DS.hess_obj[var][0] < len(DS.sVars):
                if DS.hess_obj[var][0] == 0:
                    f.write("""  // read in observations of dynamical system
  
  double number;
  fp.open("%s", ios::in | ios::binary);""" % DS.obs[var][1] + """
  if(fp.is_open()){
    while(fp >> number){
      nRead++;
      if (nRead > nSkip){
        inputData.push_back(number);
        fp.get();
      }
    }
  }
  fp.close();
  nRead = 0;
  
""")
                else: 
                    f.write("""  fp.open("%s", ios::in | ios::binary);""" % DS.obs[var][1] + """  
  if(fp.is_open()){
    while(fp >> number){
      nRead++;
      if (nRead > nSkip){
        inputData%d.push_back(number);
        fp.get();
      }
    }
  }
  fp.close();
  nRead = 0;
  
""" % ((DS.hess_obj[var][0]+1)))
                
                
        for var in range(len(DS.hess_obj)):
            if DS.hess_obj[var][0] < len(DS.sVars):
                if DS.hess_obj[var][0] == 0:
                    f.write("""  data_ = new Number[nObs_];
  for (Index i=0; i<nObs_; i++){
    data_[i] = inputData[i];
  }
""")
                else: 
                    f.write("""  data%d_ = new Number[nObs_];
  for (Index i=0; i<nObs_; i++){
    data%d_[i] = inputData%d[i];
  }
""" % ((DS.hess_obj[var][0]+1),(DS.hess_obj[var][0]+1),(DS.hess_obj[var][0]+1)))
        
        try:          
            for inj in range(len(DS.forc)):
                if inj == 0:
                    f.write("""  fp.open("%s", ios::in | ios::binary);""" % DS.forc[inj][1] + """
  if(fp.is_open()){
    while(fp >> number){
      nRead++;
      if (nRead > nSkip){
        injData.push_back(number);
        fp.get();
      }
    }
  }
  fp.close();
  nRead = 0;
  
""")
                else:
                    f.write("""  fp.open("%s", ios::in | ios::binary);""" % DS.forc[inj][1] + """
  if(fp.is_open()){
    while(fp >> number){
      nRead++;
      if (nRead > nSkip){
        injData%d.push_back(number);
        fp.get();
      }
    }
  }
  fp.close();
  nRead = 0;
  
""" % (inj+1))
                    
                    
            for inj in range(len(DS.forc)):
                if inj == 0:
                    f.write("""  inj_ = new Number[nObs_];
  for (Index i=0; i<nObs_; i++){
  inj_[i] = injData[i];
  }
""")
                else: 
                    f.write("""  inj%d_ = new Number[nObs_];
  for (Index i=0; i<nObs_; i++){
  inj%d_[i] = injData%d[i];
  }
""" % ((inj+1),(inj+1),(inj+1)))
        except:
            pass
        #the next statement ends the constructor and starts to write the destructor
        f.write("""
}
                
//destructor
%s_NLP::~%s_NLP()""" % (self.prob, self.prob) + """
{
  // make sure we delete everything we allocated
""")
                
        for var in range(len(DS.hess_obj)):
            if DS.hess_obj[var][0] < len(DS.sVars):
                if DS.hess_obj[var][0] == 0:
                    f.write("""  delete [] data_;
""")
                else: 
                    f.write("""  delete [] data%s_;
""" % (DS.hess_obj[var][0]+1))
        try:           
            for inj in range(len(DS.forc)):
                if inj == 0:
                    f.write("""  delete [] inj_;
""")
                else: 
                    f.write("""  delete [] inj%d_;
""" % (inj+1))
        except:
            pass
        #destructor: drop saved_state and bounds:
        f.write( """
  delete [] saved_state;
  delete [] bounds;
  delete [] h;
  delete [] prob;
  delete [] mssConf;
""")
                
            
        f.write("""}
        
// returns the size of the problem
bool %s_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				   Index& nnz_h_lag,
				   IndexStyleEnum& index_style)""" % self.prob + """
{
  // total number of ipopt problem variables
  n = %d*(nObs_) + %d;""" % (self.n_vars, self.n_pars) + """

  // number of problem constraints 
  m = %d*(nObs_-1)/2; """ % (self.n_constr) + """

  // number of non-zero jacobian elements
  nnz_jac_g = %d*(nObs_-1)/2;""" % (DS.nnz_jac_g) + """

  // number of non-zero elements in the hessian of the lagrangian
  nnz_h_lag = (%d + %d)*(nObs_)+%d;""" % (DS.nnz_h_con, DS.nnz_h_obj, DS.nnz_h_par) + """

  // use the C style indexing (0-based) for the matrices
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool %s_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
				      Index m, Number* g_l, Number* g_u)""" % self.prob +"""
{

  for (Index var = 0; var < n_vars; var++){
    for (Index i=var*nObs_; i<(var+1)*nObs_; i++) {
      x_l[i] = bounds[3*var];
      x_u[i] = bounds[3*var + 1];
    }
  }
		
  for (Index ctrl = 0; ctrl < n_ctrl; ctrl+=2){
    for (Index i=(n_vars+ctrl)*nObs_; i<(n_vars+ctrl+1)*nObs_; i++) {
      x_l[i] = bounds[3*(n_vars+ctrl)];
      x_u[i] = bounds[3*(n_vars+ctrl) + 1];
      x_l[i+nObs_] = bounds[3*(n_vars+ctrl+1)];
      x_u[i+nObs_] = bounds[3*(n_vars+ctrl+1) + 1];
    }
  }
		
  for (Index par = 0; par < n_pars; par++){	
    x_l[(n_vars+2*n_ctrl)*nObs_+par] = bounds[3*(n_vars+2*n_ctrl+par)];
    x_u[(n_vars+2*n_ctrl)*nObs_+par] = bounds[3*(n_vars+2*n_ctrl+par)+1];
  }        
""")
                
        f.write("""               
  // all constraints are equality constraints with right hand side zero
  for (Index j=0; j<m; j++) {
    g_l[j] = g_u[j] = 0.;
  }

  return true;
}

// returns the initial point for the problem
bool %s_NLP::get_starting_point(Index n, bool init_x, Number* x,
					 bool init_z, Number* z_L, Number* z_U,
					 Index m, bool init_lambda,
					 Number* lambda)""" % self.prob + """
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  //now using saved_state to initialise:
  cout << "Setting x from saved_state." << endl;
  for ( Index i = 0 ; i < ((n_vars+2*n_ctrl)*nObs_)+n_pars ; i++ ){
	  x[i] = saved_state[i];
  }

  //if setToData, pull x[j] to data_[j]
  if ( setToData ){
	  cout << "Setting x to data." << endl;
	  for ( Index i = 0 ; i < nObs_ ; i++ ){
		  x[i] = data_[i];
	  }
  }

  return true;
}
        
// returns the value of the objective function
bool %s_NLP::eval_f(Index n, const Number* x,
			     bool new_x, Number& obj_value)""" % self.prob + """
{
  obj_value = 0.;
  for (Index i=0; i<nObs_; i++) {
    obj_value += (%s) / objFactor;""" % self.obj_cpp + """
  }
  
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool %s_NLP::eval_grad_f(Index n, const Number* x,
				  bool new_x, Number* grad_f)""" % self.prob + """
{
  // initialise gradient of objective function
  for (Index i=0; i<n; i++){
    grad_f[i] = 0.0;
  }
  
  // specify non-zero values of grad f
  for (Index i=0; i<nObs_; i++){""") 
    
        for var in range(len(DS.grad_obj)):
            if DS.grad_obj[var][0] == 0:
                f.write("""\n    grad_f[i] = (%s) / objFactor;""" % self.convert_objective(DS, DS.grad_obj[0][1]))
            else:
                f.write("""\n    grad_f[i+%d*nObs_] = (%s) / objFactor;""" % (DS.grad_obj[var][0], self.convert_objective(DS, DS.grad_obj[var][1]))) 
            
        f.write("""
  }
        
  return true;
}
        
// return the value of the constraints: g(x)
bool %s_NLP::eval_g(Index n, const Number* x,
			     bool new_x, Index m, Number* g)""" % self.prob + """
{
  Index con = 0;
  for(Index j=0; j<(nObs_-%d); j+=2) {
  """ % self.con_limit)
        #create the array of constraint STRINGS
        #create a second array in which x[j] is replaced with data_[j]
        #set up the mss conditional
        #use array 2 on mss cases
        #else use array 1
        strings1 = []
        strings2 = []
        for var in range(len(DS.dEqns)):
            strings1.append("""     g[con+%d*(nObs_-1)/2] = %s;""" % (var,self.convert_constraint(DS, sym.ccode(DS.dEqns[var]))))
            #f.write("""\n    g[con+%d*(nObs_-1)/2] = %s;""" % (var,self.convert_constraint(DS, sym.ccode(DS.dEqns[var]))))
        for foo in strings1:
            #print( foo )
            #use data instead of x when n_mss triggers
            #and replace the zero bracket with an actual zero
            strings2.append( foo.replace( "x[j]","data_[j]" ).replace("data_[j] - data_[j]", "0.0") )
        #open the mss conditional:
        f.write("""	if ( j % n_mss == 0  and j > 0 ){
""") #have to escape the % character because format! Wow.
        #in the mss case:
        for foo in strings2:
            #print( foo )
            f.write(foo)
            f.write("""
""")
        #else:
        f.write("""	}
	else{
""")
        for foo in strings1:
            f.write(foo)
            f.write("""
""")
        #and close the conditional:
        f.write("""
    }
    con++;
  }
        
  return true;
}
        
// return the structure or values of the jacobian
bool %s_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
				 Index m, Index nele_jac, Index* iRow,
				 Index *jCol, Number* values)""" % self.prob + """
{
  if (values == NULL) {
    // return the structure of the jacobian

    Index inz = 0;
    Index con = 0;""")
        for con in range(len(DS.dEqns)):
            f.write("""
    con = 0;
    for (Index j=0; j<(nObs_-%d); j+=2) {""" % self.con_limit)
            for nnz_ele in range(len(DS.ipopt_jac_struct[con])):
                f.write("""
      iRow[inz] = con + (%s)*(nObs_-1)/2;""" % (con) + """
      jCol[inz] = %s;""" % DS.ipopt_jac_struct[con][nnz_ele] + """
      inz++;""")
            f.write("""
      con++;
    }""")
        f.write("""
    // sanity check
    assert(inz==nele_jac);

    //create & populate the correctJacobian flag array

    cout << "Finding Jacobian entries requiring MSS correction." << endl;
    //create bool array of entries that need MSS correction
    //d/dx[j] becomes zero when data_[j] replaces x[j]
    correctJacobian = new bool[nele_jac]; //array now ready
    int modulus = (nObs_-1)/2 ;//to obtain con from iRow
    int thisCon;
    for ( Index inzz = 0 ; inzz < nele_jac ; inzz++){
		correctJacobian[inzz] = false;
		if (jCol[inzz] >= nObs_ or jCol[inzz] == 0) continue; //not a d/dx[j] entry, or j=0 case
		thisCon = iRow[inzz] % modulus;
		if ( jCol[inzz] != 2*thisCon ) continue; //not d/dx[j]
		if ( jCol[inzz] % n_mss != 0 ) continue; //not an MSS point
		correctJacobian[inzz] = true; //a case for treatment
	}
    
    
  }else{
    // return the values of the jacobian of the constraints
        
    Index inz = 0;""")
        for con in range(len(DS.dEqns)):
            #handle a jacobian entry
            strings1 = []
            strings2 = []
            for nnz_ele in range(len(DS.ipopt_jac_struct[con])):
                strings1.append("""
      values[inz] = %s """ % self.convert_constraint(DS, DS.ipopt_jacobian[con][nnz_ele][1]) + """;
      inz++;
      """) #strings1 contains the unmodified equation string
            for foo in strings1:
                strings2.append( foo.replace("x[j]","data_[j]").replace("data_[j] - data_[j]","0.0") )
                #strings2 is modified for mss
            f.write("""
    for (Index j=0; j<(nObs_-%d); j+=2) {""" % self.con_limit)
            #do the mss conditional
            f.write("""	if ( j % n_mss == 0  and j > 0 ){
""")
            for foo in strings2:
                f.write( foo )
            f.write("""
       } else{
   """ )
            for foo in strings1:
                f.write( foo )
            #close both brackets:
            f.write("""
       }
""")
            f.write("""
    }""")
        #finished handling jacobian entries indexed by con
        f.write("""
    // sanity check
    assert(inz==nele_jac);

    //check the flag array and zero the correct entries
    //handle MSS corrections:
    for ( Index inzz = 0 ; inzz < nele_jac ; inzz++ ){
        if ( correctJacobian[inzz] ){
            values[inzz] = 0.0;
        }
    }
    
  }
        
  return true;
}
        
//return the structure or values of the hessian
bool %s_NLP::eval_h(Index n, const Number* x, bool new_x,
			     Number obj_factor, Index m, const Number* lambda,
			     bool new_lambda, Index nele_hess, Index* iRow,
			     Index* jCol, Number* values)""" % self.prob + """
{
  if (values == NULL) {
  
    // begin with structure of the hessian of the objective function
    Index inz = 0;
    """)    
        
        for var in range(len(DS.ipopt_hess_struct[0])):
            f.write("""for (Index j=0; j<(nObs_); j++){
      iRow[inz] = %s;""" % DS.ipopt_hess_struct[0][var][0] + """
      jCol[inz] = %s;""" % DS.ipopt_hess_struct[0][var][1] + """
      inz++;
    }
    
    """)
        f.write("""// structure of the hessian of the constraints next
    """)
        
        for con in range(1,len(DS.ipopt_hess_struct)):
            for group in range(len(DS.ipopt_hess_struct[con])): 
                if DS.ipopt_hess_struct[con][group][0][0] == 0 and len(DS.ipopt_hess_struct[con][group]) > 1:
                    f.write("""for (Index j=0; j<(nObs_); j++){
      iRow[inz] = %s;""" % DS.ipopt_hess_struct[con][group][0][2] + """
      jCol[inz] = %s;""" % DS.ipopt_hess_struct[con][group][0][3] + """
      inz++;
    """)
                    f.write("""}
    """)
                    
                    
        for con in range(1,len(DS.ipopt_hess_struct)):
            for group in range(len(DS.ipopt_hess_struct[con])): 
                if DS.ipopt_hess_struct[con][group][0][0] == 0 and len(DS.ipopt_hess_struct[con][group]) == 1:
                    f.write("""iRow[inz] = %s;""" % DS.ipopt_hess_struct[con][group][0][2] + """
            jCol[inz] = %s;""" % DS.ipopt_hess_struct[con][group][0][3] + """
            inz++;
            
        """)
        
            
        f.write("""assert(inz == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the upper right
    // triangle only
            
    // initialise hessian matrix
    Index inz = 0;
    for (Index j=0; j<(%d+%d)*(nObs_)+%d; j++){""" % (DS.nnz_h_con, DS.nnz_h_obj, DS.nnz_h_par) + """
      values[inz] = 0.0;
      inz++;
    }
            
    // begin by filling in values for the hessian of the objective function
    inz = 0;
    """)   
        obj_count = 0
        for var in range(len(DS.hess_obj)):
            f.write("""for (Index j=0; j<(nObs_); j++){
      values[inz] += obj_factor*(%s) / objFactor;""" % self.convert_constraint(DS, DS.hess_obj[var][2]) + """
      inz++;
    }                       
    """)
            obj_count += 1
                        
        f.write("""
    // next, fill in values for non-zero elements with new coordinates
    //terms involving lambda are constraint terms and so get gFactor
    Index con = 0;""")
        
        con_count = 0               
        for con in range(1,len(DS.ipopt_hess_struct)):
            for group in range(len(DS.ipopt_hess_struct[con])): 
                if DS.ipopt_hess_struct[con][group][0][0] == 0 and len(DS.ipopt_hess_struct[con][group]) > 1:
                    f.write("""
    con = 0;                
    for (Index j=0; j<(nObs_-%d); j+=2){
      inz = (%d + %d)*nObs_;""" % (self.con_limit,len(DS.ipopt_hess_struct[0]), con_count) + """
      """)
                    #introducing mss conditional:
                    strings1 = []
                    strings2 = []
                    for ele in range(len(DS.ipopt_hess_struct[con][group])):
                        strings1.append("""      values[inz+j] += lambda[con+%d*(nObs_-1)/2]*(%s);""" % (con-1,self.convert_constraint(DS, DS.ipopt_hessian[con-1][group][ele][2])) + """
      inz++;
      """)
                    for foo in strings1:
                        strings2.append( foo.replace("x[j]","data_[j]").replace("data_[j] - data_[j]","0.0") )
                    f.write("""
        if ( j % n_mss == 0  and j > 0 ){
""")
                    for foo in strings2:
                        f.write( foo )
                    f.write("""
        } else {
""")
                    for foo in strings1:
                        f.write(foo)
                    f.write("""
        }
""")
                    f.write("""con++; 
    }""")
                    con_count += 1
        
        par_count = 0
        for con in range(1,len(DS.ipopt_hess_struct)):
            for group in range(len(DS.ipopt_hess_struct[con])): 
                if DS.ipopt_hess_struct[con][group][0][0] == 0 and len(DS.ipopt_hess_struct[con][group]) == 1:
                    f.write("""
    con = 0;               
    for (Index j=0; j<(nObs_-%d); j+=2){
      """ % self.con_limit)
                    ##
                    #introducing mss conditional:
                    strings1 = []
                    strings2 = []
                    for ele in range(len(DS.ipopt_hess_struct[con][group])):
                        strings1.append("""       values[(%d+%d)*nObs_+%d] += lambda[con+%d*(nObs_-1)/2]*(%s);""" % (obj_count,con_count,par_count,con-1,self.convert_constraint(DS, DS.ipopt_hessian[con-1][group][ele][2])) + """
        """)
                    for foo in strings1:
                        strings2.append( foo.replace("x[j]","data_[j]").replace("data_[j] - data_[j]","0.0") )
                    f.write("""
        if ( j % n_mss == 0  and j > 0 ){
""")
                    for foo in strings2:
                        f.write( foo )
                    f.write("""
        } else {
""")
                    for foo in strings1:
                        f.write(foo)
                    f.write("""
        }
""")
                    ##
                    ##
                    f.write("""con++;
    }
    """)
                    par_count += 1
                
            
        f.write("""
    
    // next, fill in values for non-zero elements with used coordinates           
    """)
        #this first loop is for d/d[state] variables so this is where we look for correction:
        for con in range(1,len(DS.ipopt_hess_struct)):
            for group in range(len(DS.ipopt_hess_struct[con])): 
                if DS.ipopt_hess_struct[con][group][0][0] == 1 and len(DS.ipopt_hess_struct[con][group]) > 1:
                    f.write("""con = 0;
    for (Index j=0; j<(nObs_-%d); j+=2){
      inz = %d*nObs_;""" % (self.con_limit,DS.ipopt_hess_struct[con][group][0][1]) + """
      """)
                    ##
                    #myInt is the multiplier/offset on inz:
                    myInt = DS.ipopt_hess_struct[con][group][0][1]
                    if ( myInt == 0 ):
                        f.write( """ //in this case myInt is zero, so correct the Hessian for d/dx[j]
""" )
                    correctHessian = ( myInt == 0 ) #true in the zero case...
                    #introducing mss conditional:
                    strings1 = []
                    strings2 = []
                    for ele in range(len(DS.ipopt_hess_struct[con][group])):
                        strings1.append("""values[inz+j] += lambda[con+%d*(nObs_-1)/2]*(%s);""" % (con-1,self.convert_constraint(DS, DS.ipopt_hessian[con-1][group][ele][2])) + """
      inz++;
      """)
                        if len(DS.ipopt_hess_struct[con][group]) == 2:
                            strings1.append("""inz++;
      """)
                    for foo in strings1:
                        strings2.append( foo.replace("x[j]","data_[j]").replace("data_[j] - data_[j]","0.0") )
                    f.write("""
        if ( j % n_mss == 0  and j > 0 ){
""")
                    for foo in strings2:
                        if ( correctHessian ):
                            f.write( foo.replace( "values","//values") )
                            #simply doesn't make the change, by commenting out the values += statement
                        else:
                            f.write( foo ) #makes change to element as normal
                    f.write("""
        } else {
""")
                    for foo in strings1:
                        f.write(foo)
                    f.write("""
        }
""")
                    ##
                    ##
                    f.write("""con++;
    }     
    """)
        #this second loop is for d/d[param] so will not have any d/dx[j] correction
        for con in range(1,len(DS.ipopt_hess_struct)):
            for group in range(len(DS.ipopt_hess_struct[con])): 
                if len(DS.ipopt_hess_struct[con][group]) == 1 and DS.ipopt_hess_struct[con][group][0][0] == 1:
                    f.write("""con = 0;
    for (Index j=0; j<(nObs_-%d); j+=2){
      """ % self.con_limit)
                    ##
                    #introducing mss conditional:
                    strings1 = []
                    strings2 = []
                    for ele in range(len(DS.ipopt_hess_struct[con][group])):
                        strings1.append("""values[(%d+%d)*nObs_+%d] += lambda[con+%d*(nObs_-1)/2]*(%s);""" % (obj_count,con_count,DS.ipopt_hess_struct[con][group][0][5],
                                                con-1,self.convert_constraint(DS.ipopt_hessian[con-1][group][ele][2])) + """
      """)
                    for foo in strings1:
                        strings2.append( foo.replace("x[j]","data_[j]").replace("data_[j] - data_[j]","0.0") )
                    f.write("""
        if ( j % n_mss == 0  and j > 0 ){
""")
                    for foo in strings2:
                        f.write( foo )
                    f.write("""
        } else {
""")
                    for foo in strings1:
                        f.write(foo)
                    f.write("""
        }
""")
                    ##
                    ##
                    f.write("""con++;
    }
    """)
            
        f.write("""assert(inz == nele_hess);
  }
        
  return true;
}
        
void %s_NLP::finalize_solution(SolverReturn status,
					Index n, const Number* x,
					const Number* z_L, const Number* z_U,
					Index m, const Number* g,
					const Number* lambda,
					Number obj_value,
					const IpoptData* ip_data,
					IpoptCalculatedQuantities* ip_cq)""" % self.prob + """
{
  // here is where we store the solution to variables, parameters, write
  // to a file, etc.
  printf("With status %d", status ); //adding a report of solver status
  printf("and n_mss = %d", n_mss ); //adding a report of our current mss param
  printf(" and with objective = %.8f \\n", obj_value ); //including n_mss report
  """)
        
        for par in range(self.n_pars):
            f.write("""printf("Parameter\\t[%d]\\t[%s] """ % (par+1,DS.sPars[par]) + """ = %f\\n",""" + """ x[%d*nObs_+%d]);""" % (self.n_vars,par) + """
  """)
        f.write("""
  printf("\\nWriting objective function...\\n");
  
  FILE* objf = fopen("../outputMSS/objective_func.dat", "w");
  
  fprintf(objf, "%e", obj_value);
  
  fclose(objf);
  
  printf("\\nWriting estimated state history...\\n");
  
  FILE* svf = fopen("../outputMSS/state_history.dat", "w");
  """)
        
        sv_out_frmt = "%d\\t "
        for var in range(self.n_vars):
            sv_out_frmt += "%e\\t "
        for var in range(self.n_data):
            sv_out_frmt += "%e\\t "
        sv_out_frmt += "\\n"
        
        sv_out_vals = "i"
        for var in range(self.n_vars):
            sv_out_vals += ", x[%d*nObs_+i]" % var    
        for var in range(len(DS.hess_obj)):
            if DS.hess_obj[var][0] < len(DS.sVars):
                if DS.hess_obj[var][0] == 0:
                    sv_out_vals += ", data_[i]"
                else:
                    sv_out_vals += ", data%d_[i]" % (DS.hess_obj[var][0]+1)
        
        f.write("""
  fprintf(svf, "\\n\\nWith n_mss = %d \\n", n_mss); // reporting n_mss also...
  """)
        f.write("""
  fprintf(svf, "\\n\\nSolution of the primal variables, x\\n");
  for (Index i=0; i<nObs_; i++) {
    fprintf(svf, "%s", %s);""" % (sv_out_frmt, sv_out_vals) + """
  }
  
  fclose(svf);
  
  FILE* pf = fopen("../outputMSS/parameters.dat", "a");

  fprintf(pf, "\\n\\nWith n_mss = %d", n_mss); //including n_mss report
  fprintf(pf, " and with objective = %.8f \\n", obj_value ); //including n_mss report
  printf("\\nWriting estimated parameter values...\\n");
   
  """)
        for par in range(self.n_pars):
            f.write("""fprintf(pf, "[%d]\\t[%s]\\t""" % (par+1,DS.sPars[par]) + """%f\\t%f\\t%f\\n", """ + """bounds[3*(n_vars+2*n_ctrl+%d)], bounds[3*(n_vars+2*n_ctrl+%d)+1], x[%d*nObs_+%d])""" 
            % (par,par,self.n_vars,par) + """;
  """)
        f.write("""            
  fclose(pf);

  //add restart file saves here:
  double trackScore = obj_value / sqrt( n_mss ) ; //score scaled to favour large n_mss
  cout << "Tracking score: " << trackScore << endl;

  if ( trackScore < mssTrack ){
    //tracking file save here
    cout << "Tracking score better than previous best: " << mssTrack << endl;
    mssTrack = trackScore; //for next time
    cout << "Saving state as restart file restart.dat" << endl;
    ofstream stateFile;
    stateFile.open("../outputMSS/restart.dat");
    for ( Index j = 0 ; j < nObs_ ; j++ ){
        for ( Index i = 0 ; i < n_vars ; i++ ){
            stateFile << " " << x[ j + i*nObs_ ];
        }
        stateFile << endl;
    }
    stateFile.close();
    
    cout << "Saving restart bounds file restart.bounds.txt" << endl;
    ofstream boundsFile;
    boundsFile.open("../outputMSS/restart.bounds.txt");
    boundsFile << "#With n_mss = " << n_mss;
    boundsFile << " and obj = " << obj_value << endl;
    boundsFile << "#state variables:" << endl;
    for ( Index i = 0 ;i < n_vars ; i++ ){
        boundsFile << bounds[3*i] << "," << bounds[3*i +1] << ",";
        boundsFile << x[i*nObs_] << endl; //setting to value of first point...
    }
    boundsFile << "#control variables:" << endl;
    for ( Index i = 0 ; i < n_ctrl ; i+= 2 ){
        boundsFile << bounds[3*(n_vars + i)] << ",";
        boundsFile << bounds[3*(n_vars + i)+1] << ",";
        boundsFile << bounds[3*(n_vars + i)+2] << endl;
        boundsFile << bounds[3*(n_vars + i +1)] << ",";
        boundsFile << bounds[3*(n_vars + i +1)+1] << ",";
        boundsFile << bounds[3*(n_vars + i +1)+2]; //don't recycle control variable value
        boundsFile << endl;
    }
    boundsFile << "#parameters: " << endl;
    for ( Index i = 0 ; i < n_pars ; i++ ){
        boundsFile << bounds[3*(n_vars + 2*n_ctrl + i)] <<  ",";
        boundsFile << bounds[3*(n_vars + 2*n_ctrl + i)+1] <<  ",";
        boundsFile << x[(n_vars + 2*n_ctrl)*nObs_ + i]; //current parameter value
        boundsFile << endl;
    }
    boundsFile.close();
  }



  //using finalize_solution to update n_mss for next round:
  cout << "Current mss value: " << n_mss << endl;
  int newN = (int)( n_mss * mss_multiplier ); // multiply and cast
  //it's important that n_mss be even:
  if ( newN%2 > 0 ){
    newN += 1; //increment odd to even.
  }
  if ( newN > n_mss ){
    n_mss = newN;
  }
  else{
    n_mss += 2; //minimal increment
  }
  
  cout << "New mss value: " << n_mss << endl;
  delete [] correctJacobian; //cleanup the array!
  //save restart values from this run into saved_state

  bool passData = false; //usually using setToData anyway
  bool passCurrents = true;
  bool passControl = false; //must adapt to data model
  bool passParams = true; //required!

  if ( passData ){
    for ( Index i = 0 ; i < nObs_ ; i++ ){
      saved_state[i] = x[i];
    }
  }
  if ( passCurrents ){
    for ( Index i = nObs_ ; i < n_vars*nObs_ ; i++ ){
      saved_state[i] = x[i];
    }
  }
  if ( passControl ){
    for ( Index i = n_vars*nObs_ ; i < (n_vars+2*n_ctrl)*nObs_ ; i++ ){
      saved_state[i] = x[i];
    }
  }
  if ( passParams ){
      double verySmall = 1E-6;
    for ( Index i = (n_vars+2*n_ctrl)*nObs_ ; i < ((n_vars+2*n_ctrl)*nObs_)+n_pars ; i++ ){
      saved_state[i] = x[i];
      if (not doYank ) continue; //otherwise proceed to yank!

      //yank uses bounds information:
        //  bounds[3*(n_vars+2*n_ctrl+40)], bounds[3*(n_vars+2*n_ctrl+40)+1], x[9*nObs_+40]);
	Index whichParam = i - ((n_vars+2*n_ctrl)*nObs_);
	    double lowerBound = bounds[3*(n_vars+2*n_ctrl+whichParam)];
	    double upperBound = bounds[3*(n_vars+2*n_ctrl+whichParam)+1];
	    double range = upperBound - lowerBound;
	    if ( range < verySmall ) continue; //tiny window, pass
	if ( x[i] < lowerBound + verySmall or x[i] > upperBound - verySmall ){
	    //it's stuck at the edge of the window; reset to middle
	    saved_state[i] = 0.5*(upperBound + lowerBound);
	    cout << "YANKING parameter " << whichParam + 1;
	    cout << " from " << x[i];
	    cout << " to " << saved_state[i] << endl;
	}
      
    }
  }
  
} 
""")        
        f.close()
        
    
    def build_make(self, ipopt_path):
        """
        """
        makefile = 'Makefile'
        f = open(makefile, 'w')
        
        f.write("""# Copyright (C) 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# $Id: Makefile.in 2018

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = %s.exe""" % self.prob + """

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = %s_main.o \\""" % self.prob + """
	%s_nlp.o""" % self.prob + """

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD

# additional C++ Compiler options for linking
CXXLINKFLAGS = -Wl,--rpath -Wl,%s""" % ipopt_path + """ 

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=%s/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)""" % ipopt_path + """

# Linker flags
LIBS = `PKG_CONFIG_PATH=%s/pkgconfig: pkg-config --libs ipopt`""" % ipopt_path + """

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS) ipopt.out

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'""")
