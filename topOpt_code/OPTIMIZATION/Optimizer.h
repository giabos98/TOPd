#pragma once


#include "../CODE_HEADERS/codeHeader.h"
#include "../geometry.h"
#include "../nsProblem.h"
#include "Mean_Diffusion_Filter.h"

class OPTIMIZER
{
private:
    prec q;
    prec alpha_min;
    prec alpha_max;
    prec Vr; 
    prec V0;
    int funcType;
    int onlyGrad;
    VECTOR gamma;
    int time_integration_procedure = 0;
    //-----------
    // mesh info
    //-----------
    PHYSICS* physics;
    VECTOR_INT nodeInDom;
    VECTOR_INT elemInDom;
    VECTOR_INT optNodeFromGlobNode;
    int dim;
    int nElem_v;
    int nNodes_v;
    int n_nodes_in_dom;
    int n_elems_in_dom;
    prec rho;
    prec mu;

    std::vector<MATRIX> COEF;
    int nTimeSteps;

    prec temp_obj_functional;
    VECTOR temp_func_val;
    VECTOR temp_no_weighted_func_val;
    VECTOR temp_d_obj_functional;
    prec temp_func_in_box;
    prec temp_func_out_box;

    int first_it_flag = 1;
    prec functional_normalization_factor;

public:
    // functional weights
    VECTOR fWeights;

    // custom functional weights definition
    int customFunc;

    // functional values
    VECTOR func_val;
    VECTOR no_weighted_func_val;
    prec obj_functional;
    prec func_in_box;
    prec func_out_box;

    // functional derivative w.r.t. the optimization parameter
    VECTOR d_obj_functional;

    // inital functional value
    prec f0Init;
    VECTOR func_val_init;
    prec func_in_box_init;
    prec func_out_box_init;

    // current functional absolute value
    prec obj_abs_val;

    //alpha & dAlpha
    VECTOR alpha;
    VECTOR dAlpha;

    //smooth gamma using elemental values
    int smooth_gamma_flag = 0;

    // accelerated gamma
    int gamma_acc_case = 0;
    prec gamma_acc_mean_value = 0.5;
    prec beta_max;
    prec beta_min;
    prec beta_min_init;
    prec beta_min_threshold;
    prec beta_proj;
    int beta_interpolation;
    prec change_max;
    prec change_min;
    prec crit_change;
    prec crit_beta;
    
    prec q_beta = 0;
    VECTOR gamma_acc;
    VECTOR dgamma_acc;
    
    // diffusion filter
    int diffusion_filter_case = 0;
    Mean_DIFFUSION_FILTER diffusionFilter;
    VECTOR gamma_filter;
    VECTOR dgamma_filter;

    //constraints
    CONSTRAINTS constraints;

    //------------------------------------------------
    // CONSTRUCTOR
    //------------------------------------------------
    OPTIMIZER(){};
    //------------
    void initialize(PHYSICS* physIn, VECTOR_INT &nodeInDom_In, VECTOR_INT &elemInDom_In, VECTOR_INT &optNodeFromGlobNode_In, prec q_in, prec a_min, prec a_max, prec V0_in, prec Vr_in, int funcId, int customFunctional, int onlyGradient)
    {
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  | DEPRECATED METHOD |  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        physics = physIn;
        nodeInDom = nodeInDom_In;
        elemInDom = elemInDom_In;
        optNodeFromGlobNode = optNodeFromGlobNode_In;
        n_nodes_in_dom = nodeInDom.length;
        n_elems_in_dom = elemInDom.length;

        // int nNodes_v = (*physics).nNodes_v;
        dim = (*physics).dim;
        nElem_v = (*physics).nElem_v;
        mu = (*physics).mu;
        //----------------
        q = q_in;
        alpha_min = a_min;
        alpha_max = a_max;
        V0   = V0_in;
        Vr   = Vr_in;
        funcType = funcId;
        onlyGrad = onlyGradient;
        customFunc = customFunctional;
        //----------------
        COEF.resize(dim);
        COEF[0].define(nElem_v, dim+1, (*physics).Bloc_v.PP, (*physics).Bloc_v.P);
        COEF[1].define(nElem_v, dim+1, (*physics).Cloc_v.PP, (*physics).Cloc_v.P);
        if (dim == 3) COEF[2].define(nElem_v, dim+1, (*physics).Dloc_v.PP, (*physics).Dloc_v.P);

        func_val.setZeros(3);
        no_weighted_func_val.setZeros(3);
        temp_func_val.setZeros(3);
    }
    //---
    void initialize(PHYSICS* physIn, VECTOR_INT &nodeInDom_In, VECTOR_INT &elemInDom_In, VECTOR_INT &optNodeFromGlobNode_In, prec q_in, prec a_min, prec a_max, prec V0_in, prec Vr_in, int customFunctional, int time_integration, int onlyGradient, VECTOR beta, int opt_acceleration_case, prec betaMax, prec betaMin, int betaInterpolation, prec changeMax, prec changeMin, prec critChange, prec critBeta, int diff_filter_case, CONSTRAINTS &constraints_input, int smooth_gamma_between_element)
    {
        physics = physIn;
        nodeInDom = nodeInDom_In;
        elemInDom = elemInDom_In;
        optNodeFromGlobNode = optNodeFromGlobNode_In;
        n_nodes_in_dom = nodeInDom.length;
        n_elems_in_dom = elemInDom.length;

        // int nNodes_v = (*physics).nNodes_v;
        dim = (*physics).dim;
        nElem_v = (*physics).nElem_v;
        nNodes_v = (*physics).nNodes_v;
        rho = (*physics).rho;
        mu = (*physics).mu;
        //----------------
        q = q_in;
        alpha_min = a_min;
        alpha_max = a_max;
        V0   = V0_in;
        Vr   = Vr_in;
        time_integration_procedure = time_integration;
        fWeights = beta;
        
        define_inlet_pressure_elem();
        func_val.setZeros(fWeights.length);
        no_weighted_func_val.setZeros(fWeights.length);
        temp_func_val.setZeros(fWeights.length);
        temp_no_weighted_func_val.setZeros(fWeights.length);
        temp_d_obj_functional.setZeros(nNodes_v);
        d_obj_functional.setZeros(nNodes_v);
        alpha.setZeros(nNodes_v);
        dAlpha.setZeros(nNodes_v);
        smooth_gamma_flag = smooth_gamma_between_element;
        gamma_acc_case = opt_acceleration_case;
        gamma_acc.setZeros(n_nodes_in_dom);
        dgamma_acc.setZeros(n_nodes_in_dom);
        beta_max = betaMax;
        beta_min = betaMin;
        beta_min_init = beta_min;
        beta_min_threshold = beta_min;
        beta_proj = beta_min;
        beta_interpolation = betaInterpolation;
        change_max = changeMax;
        change_min = changeMin;
        crit_change = critChange;
        crit_beta = critBeta;
        diffusion_filter_case = diff_filter_case;
        gamma_filter.setZeros(n_nodes_in_dom);
        dgamma_filter.setZeros(n_nodes_in_dom);

        constraints = constraints_input;

        // initialize alpha and dAplha values in all the domain to their free domain values. 
        // Following this procedure only the values in the Omega_opt nodes will be updated at every optimization loop.
        for (int inod = 0; inod < nNodes_v; inod++)
        {
            alpha[inod]  = alpha_min;
            dAlpha[inod] = - (q / (q+1)) * (alpha_max - alpha_min);
        }

        customFunc = customFunctional;
        onlyGrad = onlyGradient;
        //----------------
        COEF.resize(dim);
        COEF[0].define(nElem_v, dim+1, (*physics).Bloc_v.PP, (*physics).Bloc_v.P);
        COEF[1].define(nElem_v, dim+1, (*physics).Cloc_v.PP, (*physics).Cloc_v.P);
        if (dim == 3) COEF[2].define(nElem_v, dim+1, (*physics).Dloc_v.PP, (*physics).Dloc_v.P);
    }

    ~OPTIMIZER(){};

    //------------------------------------------------
    // UPDATE NS & ADJ SOLUTIONS
    //------------------------------------------------
    void updateSol(int globIter)
    {
        nTimeSteps = globIter;
    }

    //------------------------------------------------
    // EVALUATE FUNCTIONAL AND CONSTRAINTS 
    //------------------------------------------------
    void getFunc(VECTOR x, prec &f0, MATRIX U);
    void getFuncAndDerivative(VECTOR x, prec &f0, VECTOR& df0, MATRIX U, MATRIX Ua);
    void update_val_and_derivative(VECTOR &x, prec &f0, VECTOR &df0, VECTOR &g, MATRIX& dg, prec &Vol);
    void update_val(VECTOR &x, prec &f0, VECTOR &g, prec &Vol);
    void update_volume_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr, MATRIX_INT &elem_v, VECTOR &Volume_v);
    void update_volume_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr);
    void update_subdomain_volume_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr, MATRIX_INT &elem_v, VECTOR &Volume_v);
    void update_subdomain_volume_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr);
    void update_edge_size_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr);
    void update_edge_size_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr);
    void update_discretizing_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr, MATRIX_INT &elem_v, VECTOR &Volume_v);
    void update_discretizing_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr);
    void update_WSS_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr);
    void update_WSS_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr);
    void update_constraints(VECTOR &g, MATRIX_INT &elem_v, VECTOR &Volume_v);
    void update_constraints_derivative(MATRIX &dg);
    // void decompose_solution(VECTOR &sol, MATRIX &U_sol, VECTOR &P_sol);
    // void updateVal(VECTOR &x, prec &f0, VECTOR &df0, VECTOR &g, MATRIX& dg, prec &Vol);
    // void updateJustVal(VECTOR &x, prec &f0, VECTOR &g);
    void check_gamma(VECTOR &gamma_value);
    //------------------------------------------------
    // GOC
    //------------------------------------------------
    void solveGOC(VECTOR &x, prec &Vol, prec &obj);

    //----------------
    // MMA
    //----------------
    void solveMMA(VECTOR &x, prec &Vol, prec &obj);

    //----------------
    // GCMMA
    //----------------
    void solveGCMMA(VECTOR &x, prec &Vol, prec &obj);

    prec max(prec coef1, prec coef2)
    {
        if (coef1 >= coef2) return coef1;
        else return coef2;
    }
    VECTOR max(VECTOR coef1, VECTOR coef2)
    {
        int N = coef1.length;
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] >= coef2[i]) res[i] = coef1[i];
            else res[i] = coef2[i];
        }
        return res;
    }
    VECTOR max(VECTOR coef1, prec coef2)
    {
        int N = coef1.length;
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] >= coef2) res[i] = coef1[i];
            else res[i] = coef2;
        }
        return res;
    }
    //---
    VECTOR min(prec* coef1, int coef2, int N)
    {
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] <= coef2) res[i] = coef1[i];
            else res[i] = coef2;
        }
        return res;
    }
    //---
    VECTOR min(prec* coef1, prec coef2, int N)
    {
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] <= coef2) res[i] = coef1[i];
            else res[i] = coef2;
        }
        return res;
    }
    //---
    VECTOR max(prec* coef1, int coef2, int N)
    {
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] >= coef2) res[i] = coef1[i];
            else res[i] = coef2;
        }
        return res;
    }
    //---
    VECTOR max(prec* coef1, prec coef2, int N)
    {
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] >= coef2) res[i] = coef1[i];
            else res[i] = coef2;
        }
        return res;
    }
    //---
    VECTOR max(prec* coef1, VECTOR coef2, int N)
    {
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] >= coef2[i]) res[i] = coef1[i];
            else res[i] = coef2[i];
        }
        return res;
    }
    prec min(prec coef1, prec coef2)
    {
        if (coef1 <= coef2) return coef1;
        else return coef2;
    }
    VECTOR min(VECTOR coef1, VECTOR coef2)
    {
        int N = coef1.length;
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] <= coef2[i]) res[i] = coef1[i];
            else res[i] = coef2[i];
        }
        return res;
    }
    VECTOR min(prec* coef1, VECTOR coef2, int N)
    {
        VECTOR res(N);
        for (int i = 0; i < N; i++)
        {
            if (coef1[i] <= coef2[i]) res[i] = coef1[i];
            else res[i] = coef2[i];
        }
        return res;
    }
 //-----------------------------------------
 private:
    // MMASUB
    /*
            
           This function mmasub performs one MMA-iteration, aimed at
           solving the nonlinear programming problem:
                
             Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
           subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
                       xmin_j <= x_j <= xmax_j,    j = 1,...,n
                       z >= 0,   y_i >= 0,         i = 1,...,m
        *** INPUT:
        
          m    = The number of general constraints.
          n    = The number of variables x_j.
         iter  = Current iteration number ( =1 the first time mmasub is called).
         xval  = Column vector with the current values of the variables x_j.
         xmin  = Column vector with the lower bounds for the variables x_j.
         xmax  = Column vector with the upper bounds for the variables x_j.
         xold1 = xval, one iteration ago (provided that iter>1).
         xold2 = xval, two iterations ago (provided that iter>2).
         f0val = The value of the objective function f_0 at xval.
         df0dx = Column vector with the derivatives of the objective function
                 f_0 with respect to the variables x_j, calculated at xval.
         fval  = Column vector with the values of the constraint functions f_i,
                 calculated at xval.
         dfdx  = (m x n)-matrix with the derivatives of the constraint functions
                 f_i with respect to the variables x_j, calculated at xval.
                 dfdx(i,j) = the derivative of f_i with respect to x_j.
         low   = Column vector with the lower asymptotes from the previous
                 iteration (provided that iter>1).
         upp   = Column vector with the upper asymptotes from the previous
                 iteration (provided that iter>1).
         a0    = The constants a_0 in the term a_0*z.
         a     = Column vector with the constants a_i in the terms a_i*z.
         c     = Column vector with the constants c_i in the terms c_i*y_i.
         d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
            
        *** OUTPUT:
        
         xmma  = Column vector with the optimal values of the variables x_j
                 in the current MMA subproblem.
         ymma  = Column vector with the optimal values of the variables y_i
                 in the current MMA subproblem.
         zmma  = Scalar with the optimal value of the variable z
                 in the current MMA subproblem.
         lam   = Lagrange multipliers for the m general MMA constraints.
         xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
         eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
          mu   = Lagrange multipliers for the m constraints -y_i <= 0.
         zet   = Lagrange multiplier for the single constraint -z <= 0.
          s    = Slack variables for the m general MMA constraints.
         low   = Column vector with the lower asymptotes, calculated and used
                 in the current MMA subproblem.
         upp   = Column vector with the upper asymptotes, calculated and used
                 in the current MMA subproblem.
*/
    void mmasub(VECTOR &xmma, VECTOR &ymma, prec &zmma, VECTOR &lam, VECTOR &xsi, VECTOR &eta, VECTOR &mu, prec &zet, VECTOR &s, VECTOR &low, VECTOR &upp,
                int m,int n,int iter, VECTOR xval,VECTOR xmin,VECTOR xmax,VECTOR xold1,VECTOR xold2, prec f0, 
                VECTOR df0, VECTOR g, MATRIX dg, prec a0, VECTOR a, VECTOR c,VECTOR d);



                
    /*
    This function subsolv solves the MMA subproblem:
                    
    minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
                + SUM[ ci*yi + 0.5*di*(yi)^2 ],
    
    subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
                alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
            
    Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
    Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
    */

    void subsolv(VECTOR &xmma, VECTOR &ymma, prec &zmma, VECTOR &lamma, VECTOR &xsimma, VECTOR &etamma, VECTOR &mumma, prec &zetmma, VECTOR &smma,
                int m,int n, prec epsimin, VECTOR low, VECTOR upp, VECTOR alfa,VECTOR beta, VECTOR p0,VECTOR q0,std::vector<VECTOR> P,std::vector<VECTOR> Q, 
                prec a0,VECTOR a, VECTOR b,VECTOR c,VECTOR d);

        /*    
         The left hand sides of the KKT conditions for the following
         nonlinear programming problem are calculated.
                
             Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
           subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
                       xmax_j <= x_j <= xmin_j,    j = 1,...,n
                       z >= 0,   y_i >= 0,         i = 1,...,m
        *** INPUT:
        
          m    = The number of general constraints.
          n    = The number of variables x_j.
          x    = Current values of the n variables x_j.
          y    = Current values of the m variables y_i.
          z    = Current value of the single variable z.
         lam   = Lagrange multipliers for the m general constraints.
         xsi   = Lagrange multipliers for the n constraints xmin_j - x_j <= 0.
         eta   = Lagrange multipliers for the n constraints x_j - xmax_j <= 0.
          mu   = Lagrange multipliers for the m constraints -y_i <= 0.
         zet   = Lagrange multiplier for the single constraint -z <= 0.
          s    = Slack variables for the m general constraints.
         xmin  = Lower bounds for the variables x_j.
         xmax  = Upper bounds for the variables x_j.
         df0dx = Vector with the derivatives of the objective function f_0
                 with respect to the variables x_j, calculated at x.
         fval  = Vector with the values of the constraint functions f_i,
                 calculated at x.
         dfdx  = (m x n)-matrix with the derivatives of the constraint functions
                 f_i with respect to the variables x_j, calculated at x.
                 dfdx(i,j) = the derivative of f_i with respect to x_j.
          a0   = The constants a_0 in the term a_0*z.
          a    = Vector with the constants a_i in the terms a_i*z.
          c    = Vector with the constants c_i in the terms c_i*y_i.
          d    = Vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
            
        *** OUTPUT:
        
        residu     = the residual vector for the KKT conditions.
        residunorm = sqrt(residu'*residu).
        residumax  = max(abs(residu)).
    */
    prec kktcheck(int m, int n, VECTOR x, VECTOR y, prec z, VECTOR lam, VECTOR xsi, VECTOR eta, VECTOR mu, prec zet,VECTOR s,
                 VECTOR xmin, VECTOR xmax, VECTOR df0, VECTOR g, MATRIX dg, prec a0, VECTOR a, VECTOR c, VECTOR d);

    //---------------------------------------------------------------------------------------------------------------------------
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // GCMMA METHODS
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    void asynmp(VECTOR &low, VECTOR &upp, prec &raa0, VECTOR& raa,
              int outeriter,int n,VECTOR &xval,VECTOR &xold1,VECTOR &xold2,VECTOR &xmin,VECTOR &xmax, prec raa0eps,VECTOR &raaeps,VECTOR &df0,MATRIX &dg);

    //---
    void gcmmasub(VECTOR &xmma,VECTOR &ymma,prec zmma,VECTOR &lam, VECTOR &xsi,VECTOR &eta,VECTOR &mu,prec zet, VECTOR &s, prec f0app, VECTOR &fapp,
                 int m,int n,int outeriter,prec epsimin,VECTOR &xval,VECTOR &xmin,VECTOR &xmax, VECTOR &low, VECTOR &upp, prec raa0, VECTOR &raa, prec f0,VECTOR &df0,
                 VECTOR &g,MATRIX &dg,prec a0,VECTOR &a,VECTOR & c, VECTOR &d);

    void concheck(prec &conserv, int m,prec epsimin,prec &f0app,prec &f0valnew, VECTOR &fapp, VECTOR &fvalnew);

    void raaupdate(prec &raa0, VECTOR &raa,
                   VECTOR &xmma,VECTOR &xval,VECTOR &xmin,VECTOR &xmax,VECTOR &low,VECTOR &upp,prec f0valnew,VECTOR &fvalnew,prec f0app,VECTOR &fapp, prec raa0eps,
                   VECTOR &raaeps,prec epsimin);

    //-----------------------------------------
    // DEFINE INLET ELEMENTS AND AREA
    //-----------------------------------------
    void define_inlet_pressure_elem();

public:
    void get_functional_and_opt_derivative(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U, MATRIX &Ua, MATRIX_INT &elem, VECTOR &area, VECTOR &P);

    void get_time_functional(int solution_time, MATRIX_INT &elem_v, VECTOR &Volume_v, MATRIX_INT &inlet_bound_elem, VECTOR &area_inlet_bound_elem);
    
    void eval_alpha_and_dAlpha();
    void eval_functional(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U, MATRIX_INT &elem, VECTOR &area, VECTOR &P);
    void eval_opt_derivative(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U, MATRIX &Ua);

    void eval_J_alpha(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U);
    void eval_J_gradU(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U);
    void eval_J_omega(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U);
    void eval_J_omega_by_grad(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U);
    void eval_J_p_inlet(MATRIX_INT &elem, VECTOR &area, VECTOR &P);

    void eval_gamma_acc_and_derivative();
    void eval_beta_for_projection_filter(prec temp_change);

    void topology_optimization_diffusive_filter();
};
