#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.15;

// Set the ref value for cost function
double cte_ref = 0;
double epsi_ref = 0;

//convert miles per hour to meter per second.
double v_ref = 100;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;


// state variables position label
size_t x_init = 0;
size_t y_init = x_init + N;
size_t psi_init = y_init + N;
size_t v_init = psi_init + N;
size_t cte_init = v_init + N;
size_t epsi_init = cte_init + N;

// actuator variables position label
size_t delta_init = epsi_init + N;
size_t a_init = delta_init + N - 1;



class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0;
    
    // penalize delta cte/delta epsi/delta velocity to build the cost function for further optimizing.
    for (int i = 0; i < N; i++){
        fg[0] += 5000*CppAD::pow(vars[cte_init + i] - cte_ref, 2);
        fg[0] += 2500*CppAD::pow(vars[epsi_init + i] - epsi_ref, 2);
        fg[0] += CppAD::pow(vars[v_init + i] - v_ref, 2);
    }
    // penalize the magnitude of the control input especially to steering value.
    for (int i = 0; i < N - 1; i++){
        fg[0] += 1500*CppAD::pow(vars[delta_init + i], 2);
        fg[0] += CppAD::pow(vars[a_init + i], 2);
        // add penalty for speed + steer
        fg[0] += 1100*CppAD::pow(vars[delta_init + i] * vars[v_init+i], 2);
    }
    
    // additional cost function to ensures a more smoother lane changing, especially the delta steering value.
    for (int i = 0; i < N - 2; i++){
        fg[0] += 5000 * CppAD::pow(vars[delta_init + i + 1] - vars[delta_init + i], 2);
        fg[0] += CppAD::pow(vars[a_init + i + 1] - vars[a_init + i], 2);
    }
    
    // build fg[], 1 more dimension compare to vars[]
      
    fg[x_init + 1]    = vars[x_init];
    fg[y_init + 1]    = vars[y_init];
    fg[psi_init + 1]  = vars[psi_init];
    fg[v_init + 1]    = vars[v_init];
    fg[cte_init + 1]  = vars[cte_init];
    fg[epsi_init + 1] = vars[epsi_init];
      
    // Define constrains
      
    for (int i = 0; i < N - 1; i++){
        // The state at time t+1 .
        AD<double> x1    = vars[x_init + i + 1];
        AD<double> y1    = vars[y_init + i + 1];
        AD<double> psi1  = vars[psi_init + i + 1];
        AD<double> v1    = vars[v_init + i + 1];
        AD<double> cte1  = vars[cte_init + i + 1];
        AD<double> epsi1 = vars[epsi_init + i + 1];
        
        // The state at time t.
        AD<double> x0    = vars[x_init + i];
        AD<double> y0    = vars[y_init + i];
        AD<double> psi0  = vars[psi_init + i];
        AD<double> v0    = vars[v_init + i];
        AD<double> cte0  = vars[cte_init + i];
        AD<double> epsi0 = vars[epsi_init + i];
        
        // Only consider the actuation at time t.
        AD<double> delta0 = vars[delta_init + i];
        AD<double> a0     = vars[a_init + i];
        
        
        // prepare for the motion model
        AD<double> f0      = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0,2) + coeffs[3] * CppAD::pow(x0,3);
        AD<double> psides0 = CppAD::atan(coeffs[1] + (2 * coeffs[2] * x0) + (3 * coeffs[3]* CppAD::pow(x0,2)));
        
        // rewrite motion model into constrains and save to fg{}.
        fg[2 + x_init + i]        = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
        fg[2 + y_init + i]        = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
        fg[2 + psi_init + i]      = psi1 - (psi0 + v0 * delta0 / Lf * dt);
        fg[2 + v_init + i]        = v1 - (v0 + a0 * dt);
        fg[2 + cte_init + i]      = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
        fg[2 + epsi_init + i]     = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
        
      }
    
  
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;

  typedef CPPAD_TESTVECTOR(double) Dvector;
    
  // state vector from latency stransformation
  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];
    
    
  

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  //size_t n_vars = 0;
    
  // number of independent variables
  // N timestep = (N - 1) timeslots
  size_t n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  // number of constraints
  // size_t n_constraints = 0;
  size_t n_constraints = N * 6;
    
  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
    
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
   
    
  // set initial variable values
  vars[x_init] = x;
  vars[y_init] = y;
  vars[psi_init] = psi;
  vars[v_init] = v;
  vars[cte_init] = cte;
  vars[epsi_init] = epsi;
    
  // set lower and upper limit
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

    
  // Genearlize all the variable bounds
  for(int i=0;i<delta_init;i++){
      vars_lowerbound[i] = -1.0e19;
      vars_upperbound[i] = 1.0e19;
  }
  
  // The upper and lower limits of delta are set to -25 and 25 degrees (values in radians).
  for (int i = delta_init; i < a_init; i++) {
      vars_lowerbound[i] = -0.436332;
      vars_upperbound[i] = 0.436332;
  }
    
    
  // Acceleration/decceleration upper and lower limits [-1.0,1.0].
  for (int i = a_init; i < n_vars; i++) {
      vars_lowerbound[i] = -1.0;
      vars_upperbound[i] =  1.0;
 }
    

  // TODO: Set lower and upper limits for variables.

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  //Lower and upper limits for constraints, All of these should be 0 except the initial state indices.
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
    
  constraints_lowerbound[x_init] = x;
  constraints_lowerbound[y_init] = y;
  constraints_lowerbound[psi_init] = psi;
  constraints_lowerbound[v_init] = v;
  constraints_lowerbound[cte_init] = cte;
  constraints_lowerbound[epsi_init] = epsi;
    
  constraints_upperbound[x_init] = x;
  constraints_upperbound[y_init] = y;
  constraints_upperbound[psi_init] = psi;
  constraints_upperbound[v_init] = v;
  constraints_upperbound[cte_init] = cte;
  constraints_upperbound[epsi_init] = epsi;
    
  

  // object that computes objective and constraints
//  vector<double> previous_actuations = {delta_prev,a_prev};
    
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
//  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
    
    
  // push back the model forecasted trajectory(green path)
  mpc_x_vals = {};
  mpc_y_vals = {};
  for (int i = 0; i < N - 1; i++) {
      mpc_x_vals.push_back(solution.x[x_init + i + 1]);
      mpc_y_vals.push_back(solution.x[y_init + i + 1]);
//      std::cout << "mpc_x_vals" <<  solution.x[x_init + i + 1] << std::endl;
//      std::cout << "mpc_y_vals " << solution.x[y_init + i + 1] << std::endl;
  }
  // return optimized actuator
  return {solution.x[delta_init], solution.x[a_init]};
}
