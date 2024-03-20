#include "OPTIMIZER.h"

void OPTIMIZER::getFunc(VECTOR x, prec &f0, MATRIX U)
{
    f0 = 0;
    int nElem_v = (*physics).nElem_v;
    MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
    VECTOR Volume_v = (*physics).Volume_v;

    prec factor = q*(alpha_min - alpha_max) + alpha_min;
    VECTOR alpha = (x-1).pointdiv(x+q);
    alpha *= factor; 
    switch (funcType)
    {
        case 0: // func = \alpha ||u||^2
        {
            for (int iel = 0; iel < n_elems_in_dom; iel++)
            {
                int globEl = elemInDom[iel];
                for (int iloc = 0; iloc < dim+1; iloc++)
                {
                    int iglob = elem_v[globEl][iloc];
                    int optNode = optNodeFromGlobNode[iglob];
                    prec tempFactor = 0;
                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        prec Uval = U[icomp][iglob];
                        tempFactor += Uval*Uval;
                    }
                    f0 += alpha[optNode] * tempFactor/(dim+1)*Volume_v[globEl]; 
                }
            }
            break;
        }
        case 1: //func = \alpha ||u||^2 + 1/2*mu ||grad(u) + grad(u)^T||^2
        {
            for (int iel = 0; iel < n_elems_in_dom; iel++)
            {
                int globEl = elemInDom[iel];
                for (int iloc = 0; iloc < dim+1; iloc++)
                {
                    int iglob = elem_v[globEl][iloc];
                    int optNode = optNodeFromGlobNode[iglob];
                    prec tempFactor = 0;
                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        prec Uval = U[icomp][iglob];
                        tempFactor += Uval*Uval;
                    }
                    f0 += alpha[optNode] * tempFactor/(dim+1)*Volume_v[globEl]; 
                }
            }
            prec part2 = 0;
            for (int iel = 0; iel < n_elems_in_dom; iel++)
            {
                MATRIX tempGrad;
                tempGrad.zeros(dim,dim);
                int globEl = elemInDom[iel];
                for (int iloc = 0; iloc < dim+1; iloc++)
                {
                    int iglob = elem_v[globEl][iloc];
                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        for (int jcomp = 0; jcomp < dim; jcomp++)
                        {
                            tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[icomp][globEl][iloc] + U[jcomp][iglob]*COEF[jcomp][globEl][iloc];
                        }
                    }
                }

                prec tempNorm = tempGrad.normFro();
                part2 += tempNorm*tempNorm * Volume_v[globEl];
            }
            part2 *= 0.5 * mu;
            f0 += part2;
            break;
        }
        case 2: // 1/2*mu ||grad(u) + grad(u)^T||^2
        {
            for (int iel = 0; iel < n_elems_in_dom; iel++)
            {
                int globEl = elemInDom[iel];
                MATRIX tempGrad;
                tempGrad.zeros(dim,dim);
                for (int iloc = 0; iloc < dim+1; iloc++)
                {
                    int iglob = elem_v[globEl][iloc];

                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        for (int jcomp = 0; jcomp < dim; jcomp++)
                        {
                            tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[icomp][globEl][iloc] + U[jcomp][iglob]*COEF[jcomp][globEl][iloc];
                        }
                    }
                }

                prec tempNorm = tempGrad.normFro();
                f0 += tempNorm*tempNorm * Volume_v[globEl];
            }
            f0 *= 0.5 * mu;
            break;
        }
    }
}

void OPTIMIZER::getFuncAndDerivative(VECTOR x, prec &f0, VECTOR& df0, MATRIX U, MATRIX Ua)
{
    f0 = 0;
    df0.resetZeros();
    int nElem_v = (*physics).nElem_v;
    MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
    VECTOR Volume_v = (*physics).Volume_v;
    prec rho = (*physics).rho;

    prec factor = q*(alpha_min - alpha_max) + alpha_min;
    VECTOR alpha = (x-1).pointdiv(x+q);
    alpha *= factor;
    if (customFunc == 0) 
    {
        switch (funcType)
        {
            case 0: // func = \alpha ||u||^2
            {
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    int globEl = elemInDom[iel];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];
                        int optNode = optNodeFromGlobNode[iglob];
                        prec tempFactor = 0;
                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            prec Uval = U[icomp][iglob];
                            tempFactor += Uval*Uval;
                        }
                        f0 += alpha[optNode] * tempFactor/(dim+1)*Volume_v[globEl]; 
                    }
                }
                break;
            }
            case 1: //func = \alpha ||u||^2 + 1/2*mu ||grad(u) + grad(u)^T||^2
            {
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    int globEl = elemInDom[iel];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];
                        int optNode = optNodeFromGlobNode[iglob];
                        prec tempFactor = 0;
                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            prec Uval = U[icomp][iglob];
                            tempFactor += Uval*Uval;
                        }
                        f0 += alpha[optNode] * tempFactor/(dim+1)*Volume_v[globEl]; 
                    }
                }
                prec part2 = 0;
                prec gradFact = 1;
                if (onlyGrad == 1) gradFact = 0;
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    MATRIX tempGrad;
                    tempGrad.zeros(dim,dim);
                    int globEl = elemInDom[iel];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];
                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            for (int jcomp = 0; jcomp < dim; jcomp++)
                            {
                                tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[icomp][globEl][iloc] + gradFact*U[jcomp][iglob]*COEF[jcomp][globEl][iloc];
                            }
                        }
                    }

                    prec tempNorm = tempGrad.normFro();
                    part2 += tempNorm*tempNorm * Volume_v[globEl];
                }
                part2 *= 0.5 * mu;
                f0 += part2;
                break;
            }
            case 2: // 1/2*mu ||grad(u) + grad(u)^T||^2
            {
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    int globEl = elemInDom[iel];
                    MATRIX tempGrad;
                    tempGrad.zeros(dim,dim);
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];

                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            for (int jcomp = 0; jcomp < dim; jcomp++)
                            {
                                tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[icomp][globEl][iloc] + U[jcomp][iglob]*COEF[jcomp][globEl][iloc];
                            }
                        }
                    }

                    prec tempNorm = tempGrad.normFro();
                    f0 += tempNorm*tempNorm * Volume_v[globEl];
                }
                f0 *= 0.5 * mu;
                break;
            }
        }
    }
    else if (customFunc == 1)
    {
        // alpha
        prec part0 = 0;
        for (int iel = 0; iel < n_elems_in_dom; iel++)
        {
            int globEl = elemInDom[iel];
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[globEl][iloc];
                int optNode = optNodeFromGlobNode[iglob];
                prec tempFactor = 0;
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    prec Uval = U[icomp][iglob];
                    tempFactor += Uval*Uval;
                }
                part0 += alpha[optNode] * tempFactor/(dim+1)*Volume_v[globEl]; 
            }
        }
        // std::cout << "\npeso: " << fWeights[0] << "\n";
        // std::cout << "\n alpha*||u||^2: " << part0 << " ";
        part0 *= fWeights[0];
        temp_func_val[0] = part0;
        // std::cout << "\n beta(1)*alpha*||u||^2: " << part0 << " ";
        // grad
        prec part1 = 0;
        prec gradFact = 1;
        if (onlyGrad == 1) gradFact = 0;
        for (int iel = 0; iel < n_elems_in_dom; iel++)
        {
            MATRIX tempGrad;
            tempGrad.zeros(dim,dim);
            int globEl = elemInDom[iel];
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[globEl][iloc];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    for (int jcomp = 0; jcomp < dim; jcomp++)
                    {
                        tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[icomp][globEl][iloc] + gradFact*U[jcomp][iglob]*COEF[jcomp][globEl][iloc];
                    }
                }
            }
            prec tempNorm = tempGrad.normFro();
            part1 += tempNorm*tempNorm * Volume_v[globEl];
        }
        part1 *= 0.5*mu;
        // std::cout << "\npeso1: " << fWeights[1] << "\n";
        // std::cout << "\n1/2*mu*||grad(u) + grad(u)^T||^2: " << part1 << " ";
        part1 *= fWeights[1];
        temp_func_val[1] = part1;
        // std::cout << "\n beta(2)*1/2*mu*||grad(u) + grad(u)^T||^2: " << part1 << " ";
        // update f0
        f0 += part0 + part1;
        //curl
        prec part2 = 0;
        VECTOR_INT ezComps;
        if (dim == 3)
        {
            ezComps.initialize(5);
            ezComps[0] = 0; ezComps[1] = 1; ezComps[2] = 2; ezComps[3] = 0; ezComps[4] = 1;
        }

        for (int iel = 0; iel < n_elems_in_dom; iel++)
        {
            VECTOR tempCurl(3);
            tempCurl.resetZeros();
            int globEl = elemInDom[iel];
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[globEl][iloc];
                if ( dim == 2)
                {
                    tempCurl[0] = 0; tempCurl[1] = 0;
                    tempCurl[2] += U[1][iglob] * COEF[0][globEl][iloc] - U[0][iglob]*COEF[1][globEl][iloc];
                }
                else if (dim == 3)
                {
                    for (int icomp = 0; icomp < dim; icomp++)
                    {
                        int comp2 = ezComps[icomp+2];
                        int comp1 = ezComps[icomp+1];
                        tempCurl[icomp] += U[comp2][iglob] * COEF[comp1][globEl][iloc] - U[comp1][iglob]*COEF[comp2][globEl][iloc];
                    }
                }
            }
            prec tempNorm = tempCurl.norm();
            part2 += tempNorm*tempNorm * Volume_v[globEl]; 
        }
        part2 *= 0.5*mu;
        // std::cout << "\npeso1: " << fWeights[2] << "\n";
        // std::cout << "\n1/2*mu*||curl(u)||^2: " << part2 << " ";
        part2 *= fWeights[2];
        temp_func_val[2] = part2;
        // std::cout << "\n beta(3)*1/2*mu*||curl(u)||^2: " << part2 << "\n";
        // update f0
        f0 += part2;
    }
    
    //--------------------------------------------------------------
    // DERIVATIVE OF THE FUCNTIONAL df0/dx
    VECTOR dAlpha(n_nodes_in_dom);
    for (int iglob = 0; iglob < n_nodes_in_dom; iglob++)
    {
        dAlpha[iglob] = -(alpha_max - alpha_min) * q * (q+1) / ((q+x[iglob])*(q+x[iglob]));
    }
    if (customFunc == 0)
    {
        switch (funcType)
        {
            case 1: // func = \alpha ||u||^2
            {
                /* PRIMA PARTE
                    d_alpha/d_x = - (q+1)/(q+x)^2
                    1)--> d_A / d_x   = d_alpha/d_x * \phi_i * ||u||^2
                %%% SECONDA PARTE
                    2)--> d_alpha/d_x phi_i u * u_a

                %%% TOTALE
                    % 1) + 2) = (d_alpha/d_x phi_i u) * (u + u_a)
                    */

                VECTOR tempU(dim);
                VECTOR tempUa(dim);
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    int globEl = elemInDom[iel];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];
                        int optNode = optNodeFromGlobNode[iglob];
                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            tempU[icomp] = U[icomp][iglob];
                            tempUa[icomp] = Ua[icomp][iglob];
                        }
                        VECTOR tempSum = tempU + tempUa*(1/rho);
                        df0[optNode] += (dAlpha[optNode] * tempU.dot(tempSum)) / (dim+1) * Volume_v[globEl];
                    }
                }
                break;
            }
            case 2:
            {
                VECTOR tempU(dim);
                VECTOR tempUa(dim);
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    int globEl = elemInDom[iel];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];
                        int optNode = optNodeFromGlobNode[iglob];
                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            tempU[icomp] = U[icomp][iglob];
                            tempUa[icomp] = Ua[icomp][iglob];
                        }
                        VECTOR tempSum = tempU + tempUa*(1/rho);
                        df0[optNode] += (dAlpha[optNode] * tempU.dot(tempSum)) / (dim+1) * Volume_v[globEl];
                    }
                }
                break;
            }
            case 3: // func = mu ||grad(u) + grad(u)^T||^2
            {
                VECTOR tempU(dim);
                VECTOR tempUa(dim);
                for (int iel = 0; iel < n_elems_in_dom; iel++)
                {
                    int globEl = elemInDom[iel];
                    for (int iloc = 0; iloc < dim+1; iloc++)
                    {
                        int iglob = elem_v[globEl][iloc];
                        int optNode = optNodeFromGlobNode[iglob];
                        for (int icomp = 0; icomp < dim; icomp++)
                        {
                            tempU[icomp] = U[icomp][iglob];
                            tempUa[icomp] = Ua[icomp][iglob];
                        }

                        df0[optNode] += (dAlpha[optNode] * tempU.dot(tempUa)) / (dim+1) * Volume_v[globEl];
                    }
                }
                break;
            }
        }
    }
    else if (customFunc == 1)
    {
        VECTOR tempU(dim);
        VECTOR tempUa(dim);
        
        for (int iel = 0; iel < n_elems_in_dom; iel++)
        {
            int globEl = elemInDom[iel];
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[globEl][iloc];
                int optNode = optNodeFromGlobNode[iglob];
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    tempU[icomp] = U[icomp][iglob];
                    tempUa[icomp] = Ua[icomp][iglob];
                }
                df0[optNode] += 1/rho*(dAlpha[optNode] * tempU.dot(tempUa)) / (dim+1) * Volume_v[globEl];
                df0[optNode] += fWeights[0]*(dAlpha[optNode] * tempU.dot(tempU)) / (dim+1) * Volume_v[globEl];
            }
        }
    }
    else throw_line("\nERROR: customFunc != 0 or 1\n");
}

// void OPTIMIZER::updateJustVal(VECTOR &x, prec &f0,  VECTOR &g)
// {
//     f0 = 0;
//     int nElem_v = (*physics).nElem_v;
//     int nNode_v = (*physics).nNodes_v;

//     MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
//     VECTOR Volume_v = (*physics).Volume_v;

//     prec factor = q*(alpha_min - alpha_max) + alpha_min;
//     VECTOR alpha = (x-1).pointdiv(x+q);
//     alpha *= factor; 
    
//     //--------------------------------------------------------------
//     // f_0, FUNCTIONAL DEFINITION
//     //--------------------------------------------------------------
//     MATRIX U(dim, nNode_v);
//     if ( (*physics).isStationary == 1) 
//     {
//         std::string currNSName = "NSCurrSol/0.txt";
//         std::string currADJName = "ADJCurrSol/0.txt";

//         std::ifstream rf;
//         rf.open(&currNSName[0], std::ios::out | std::ios::binary);
//         if(!rf) {
//             std::cout << "Cannot open file!" << std::endl;
//             throw_line("ERROR: Can't open NS sol file\n");
//         }
//         prec deltaT;
//         rf.read((char *) &deltaT, sizeof(prec));

//         for (int icomp = 0; icomp < dim; icomp++) 
//         { for (int i = 0; i < nNode_v; i++) rf.read((char *) &U[icomp][i], sizeof(prec));}
//         rf.close();
//         //---------
//         rf.open(&currADJName[0], std::ios::out | std::ios::binary);
//         if(!rf) {
//             std::cout << "Cannot open file!" << std::endl;
//             throw_line("ERROR: Can't open ADJ sol file\n");
//         }
//         rf.close();
//         getFunc(x, f0, U);
//     } 
//     else
//     {
//         prec totTime = (*physics).t_end;

//         for (int itime = 0; itime < nTimeSteps+1; itime++)
//         {
//             std::string currNSName = "NSCurrSol/" + std::to_string(itime) + ".txt";
//             std::string currADJName = "ADJCurrSol/" + std::to_string(itime) + ".txt";
//             std::ifstream rf;
//             rf.open(&currNSName[0], std::ios::out | std::ios::binary);
//             if(!rf) {
//                 std::cout << "Cannot open file!" << std::endl;
//                 throw_line("ERROR: Can't open NS sol file\n");
//             }
//             prec deltaT;
//             rf.read((char *) &deltaT, sizeof(prec));

//             MATRIX U(dim, nNode_v);
//             for (int icomp = 0; icomp < dim; icomp++) 
//             { for (int i = 0; i < nNode_v; i++) rf.read((char *) &U[icomp][i], sizeof(prec));}
//             rf.close();
//             //---------
//             rf.open(&currADJName[0], std::ios::out | std::ios::binary);
//             if(!rf) {
//                 std::cout << "Cannot open file!" << std::endl;
//                 throw_line("ERROR: Can't open ADJ sol file\n");
//             }
//             MATRIX Ua(dim, nNode_v);
//             for (int icomp = 0; icomp < dim; icomp++) 
//             { for (int i = 0; i < nNode_v; i++) rf.read((char *) &Ua[icomp][i], sizeof(prec));}
//             rf.close();
//             prec tempF0;
//             getFunc(x, tempF0, U);
//             f0 += tempF0*deltaT/2; 
//         }
//         f0 /= totTime; 
//     }

//     if (first_it_flag == 1)
//     {
//         f0Init = f0;
//         functional_normalization_factor = abs(f0Init);
//         first_it_flag = 0;
//     }
//     (*physics).f0Init = f0Init;
//     obj_abs_val = f0;
//     f0 /= functional_normalization_factor; 
//     // // --------------------------------------------------------------
//     // // MESH INDIPENDENCY FILTER (MAYBE MOVE IT AWAY PLEASE)
//     // // --------------------------------------------------------------
//     // VECTOR newDf0(n_nodes_in_dom); newDf0.resetZeros();
//     // VECTOR countW(n_nodes_in_dom); countW.resetZeros();
//     // for (int iel = 0; iel < n_elems_in_dom; iel++)
//     // {
//     //     int globEl = elemInDom[iel];
//     //      for (int iloc = 0; iloc < dim+1; iloc++)
//     //     {
//     //         int iglob = elem_v[globEl][iloc];
//     //         int optNode = optNodeFromGlobNode[iglob];
//     //         newDf0[optNode] += df0[optNode]*Volume_v[globEl]/(dim+1);
//     //         countW[optNode] += Volume_v[globEl]/(dim+1);
//     //     }
//     // }
//     // newDf0 /= countW;
//     // df0 = newDf0;

//     // newDf0.resetZeros();
//     // countW.resetZeros();
//     // for (int iel = 0; iel < n_elems_in_dom; iel++)
//     // {
//     //     int globEl = elemInDom[iel];
//     //      for (int iloc = 0; iloc < dim+1; iloc++)
//     //     {
//     //         int iglob = elem_v[globEl][iloc];
//     //         int optNode = optNodeFromGlobNode[iglob];
//     //         newDf0[optNode] += df0[optNode]*Volume_v[globEl]/(dim+1);
//     //         countW[optNode] += Volume_v[globEl]/(dim+1);
//     //     }
//     // }
//     // newDf0 /= countW;
//     // df0 = newDf0;
//     //--------------------------------------------------------------
//     // CONSTRAINTS DEFINITION
//     //--------------------------------------------------------------
//     // volume constraint
//     //--------------------
//     // int_\Omega { x } - V_r*V_0 <= 0 
//     prec integral = 0;
//     for (int iel = 0; iel < n_elems_in_dom; iel++)
//     {
//         int globEl = elemInDom[iel];
//         for (int iloc = 0; iloc < dim+1; iloc++)
//         {
//             int iglob = elem_v[globEl][iloc];
//             int optNode = optNodeFromGlobNode[iglob];
//             integral += x[optNode]/(dim+1)*Volume_v[globEl];
//         }
//     }
//     g[0]  = (integral - Vr * V0)/V0;
// }

// void OPTIMIZER::updateVal(VECTOR x, prec &f0, VECTOR &df0, VECTOR &g, MATRIX& dg, prec &Vol)
// {
//     obj_functional = 0.0;
//     func_in_box = 0.0;
//     f0 = 0.0;
    
//     func_val.resetZeros();

//     d_obj_functional.resetZeros();
//     df0.resetZeros();
//     int nElem_v = (*physics).nElem_v;
//     int nNode_v = (*physics).nNodes_v;
//     int n_node   = (*physics).nNodes;
//     int n_in_bd_elem = (*physics).n_inlet_bounds_elems;

//     MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
//     VECTOR Volume_v = (*physics).Volume_v;
//     MATRIX_INT inlet_bound_elem(n_in_bd_elem, dim, (*physics).inlet_bounds_elems.PP, (*physics).inlet_bounds_elems.P);
//     VECTOR area_inlet_bound_elem((*physics).area_inlet_bounds_elems.P, n_in_bd_elem);

//     // prec factor = q*(alpha_min - alpha_max) + alpha_min;
//     // VECTOR alpha = (x-1).pointdiv(x+q);
//     // alpha *= factor; 
    
//     //--------------------------------------------------------------
//     // f_0, FUNCTIONAL DEFINITION
//     //--------------------------------------------------------------
//     MATRIX U(dim, nNode_v);
//     VECTOR P(n_node);
//     MATRIX Ua(dim, nNode_v);
//     VECTOR Pa(n_node);

//     if ( (*physics).isStationary == 1) 
//     {
//         std::string currNSName = "NSCurrSol/0.txt";
//         std::string currADJName = "ADJCurrSol/0.txt";

//         std::ifstream rf;
//         rf.open(&currNSName[0], std::ios::out | std::ios::binary);
//         if(!rf) {
//             std::cout << "Cannot open file!" << std::endl;
//             throw_line("ERROR: Can't open NS sol file\n");
//         }
//         prec deltaT;
//         rf.read((char *) &deltaT, sizeof(prec));
//         for (int icomp = 0; icomp < dim; icomp++) 
//         { 
//             for (int i = 0; i < nNode_v; i++)
//             {
//                 rf.read((char *) &U[icomp][i], sizeof(prec));
//             } 
//         }
//         for (int inode = 0; inode < n_node; inode++)
//         {
//             rf.read((char *) &P[inode], sizeof(prec));
//         }
//         rf.close();
//         //---------
//         rf.open(&currADJName[0], std::ios::out | std::ios::binary);
//         if(!rf) {
//             std::cout << "Cannot open file!" << std::endl;
//             throw_line("ERROR: Can't open ADJ sol file\n");
//         }
//         for (int icomp = 0; icomp < dim; icomp++) 
//         { 
//             for (int i = 0; i < nNode_v; i++)
//             {
//                 rf.read((char *) &Ua[icomp][i], sizeof(prec));
//             }
//         }
//         for (int inode = 0; inode < n_node; inode++)
//         {
//             rf.read((char *) &Pa[inode], sizeof(prec));
//         }
//         rf.close();
//         temp_func_val.resetZeros();
//         getFuncAndDerivative(x, f0, df0, U, Ua);
//         func_val = temp_func_val;
//     } 
//     else
//     {
//         prec totTime = (*physics).t_end;

//         for (int itime = 0; itime < nTimeSteps+1; itime++)
//         {
//             std::string currNSName = "NSCurrSol/" + std::to_string(itime) + ".txt";
//             std::string currADJName = "ADJCurrSol/" + std::to_string(itime) + ".txt";
//             std::ifstream rf;

//             //--------
//             // read NS
//             rf.open(&currNSName[0], std::ios::out | std::ios::binary);
//             if(!rf) {
//                 std::cout << "Cannot open file!" << std::endl;
//                 throw_line("ERROR: Can't open NS sol file\n");
//             }
//             prec deltaT;
//             rf.read((char *) &deltaT, sizeof(prec));

//             MATRIX U(dim, nNode_v);
//             VECTOR P(n_node);
//             for (int icomp = 0; icomp < dim; icomp++) 
//             { 
//                 for (int i = 0; i < nNode_v; i++)
//                 { 
//                     rf.read((char *) &U[icomp][i], sizeof(prec));
//                 }
//             }
//             for (int inode = 0; inode < n_node; inode++)
//             {
//                 rf.read((char *) &P[inode], sizeof(prec));
//             }
//             rf.close();

//             //--------
//             // read ADJ
//             rf.open(&currADJName[0], std::ios::out | std::ios::binary);
//             if(!rf) {
//                 std::cout << "Cannot open file!" << std::endl;
//                 throw_line("ERROR: Can't open ADJ sol file\n");
//             }

//             MATRIX Ua(dim, nNode_v);
//             VECTOR Pa(n_node);
//             for (int icomp = 0; icomp < dim; icomp++) 
//             { 
//                 for (int i = 0; i < nNode_v; i++)
//                 {
//                    rf.read((char *) &Ua[icomp][i], sizeof(prec)); 
//                 } 
//             }
//             for (int inode = 0; inode < n_node; inode++)
//             {
//                 rf.read((char *) &Pa[inode], sizeof(prec));
//             }
//             rf.close();
//             temp_obj_functional = 0.0;
//             temp_func_val.resetZeros();
//             temp_d_obj_functional.resetZeros();
//             // prec tempF0;
//             // VECTOR tempDF0(n_nodes_in_dom);
//             //getFuncAndDerivative(x, tempF0, tempDF0, U, Ua);
//             get_functional_and_opt_derivative(x, elem_v, Volume_v, U, Ua, inlet_bound_elem, area_inlet_bound_elem, P);
//             // std::cout << temp_obj_functional << "\n";
//             // pause();
//             obj_functional += temp_obj_functional*deltaT/2; 
//             func_in_box += temp_func_in_box*deltaT/2;
//             d_obj_functional += temp_d_obj_functional*deltaT/2;
//             func_val += (temp_func_val * (deltaT / 2));
//         }
//         obj_functional /= totTime; 
//         func_in_box /= totTime;
//         d_obj_functional /= totTime;
//         func_val = func_val / totTime;
//     }

    
//     if (first_it_flag == 1)
//     {
//         f0Init = obj_functional;
//         func_val_init = func_val;
//         func_in_box_init = func_in_box;
//         func_out_box_init = f0Init - func_in_box_init;
//         if (f0Init < 1e-10)
//         {
//             functional_normalization_factor = 1;
//         } 
//         else 
//         {
//             functional_normalization_factor = abs(f0Init);
//         }
//         (*physics).f0Init = f0Init;
//         first_it_flag = 0;
//     }
    
//     obj_abs_val = obj_functional;
//     obj_functional /= functional_normalization_factor; 
//     func_in_box /= functional_normalization_factor;
//     func_out_box = obj_functional - func_in_box;
//     func_val = func_val / functional_normalization_factor;
//     d_obj_functional /= functional_normalization_factor;
//     // std::cout << "box: "<< func_in_box << "\n";
//     // std::cout << "tot: "<< obj_functional << "\n";

//     // alpha.printRowMatlabInLine("alpha");
//     // dAlpha.printRowMatlabInLine("dalpha");
//     // d_obj_functional.printRowMatlabInLine("df");
//     // std::cout << "func: " << obj_functional << "\n";
//     // func_val.printRowMatlabInLine("func val");
//     //pause();

//     f0 = obj_functional;
//     for (int ioptnode = 0; ioptnode < n_nodes_in_dom; ioptnode++)
//     {
//         int iglob = nodeInDom[ioptnode];
//         df0[ioptnode] = d_obj_functional[iglob];
//     }
//     // // --------------------------------------------------------------
//     // // MESH INDIPENDENCY FILTER (MAYBE MOVE IT AWAY PLEASE)
//     // // --------------------------------------------------------------
//     // VECTOR newDf0(n_nodes_in_dom); newDf0.resetZeros();
//     // VECTOR countW(n_nodes_in_dom); countW.resetZeros();
//     // for (int iel = 0; iel < n_elems_in_dom; iel++)
//     // {
//     //     int globEl = elemInDom[iel];
//     //      for (int iloc = 0; iloc < dim+1; iloc++)
//     //     {
//     //         int iglob = elem_v[globEl][iloc];
//     //         int optNode = optNodeFromGlobNode[iglob];
//     //         newDf0[optNode] += df0[optNode]*Volume_v[globEl]/(dim+1);
//     //         countW[optNode] += Volume_v[globEl]/(dim+1);
//     //     }
//     // }
//     // newDf0 /= countW;
//     // df0 = newDf0;

//     // newDf0.resetZeros();
//     // countW.resetZeros();
//     // for (int iel = 0; iel < n_elems_in_dom; iel++)
//     // {
//     //     int globEl = elemInDom[iel];
//     //      for (int iloc = 0; iloc < dim+1; iloc++)
//     //     {
//     //         int iglob = elem_v[globEl][iloc];
//     //         int optNode = optNodeFromGlobNode[iglob];
//     //         newDf0[optNode] += df0[optNode]*Volume_v[globEl]/(dim+1);
//     //         countW[optNode] += Volume_v[globEl]/(dim+1);
//     //     }
//     // }
//     // newDf0 /= countW;
//     // df0 = newDf0;

//     // newDf0.resetZeros();
//     // countW.resetZeros();
//     // for (int iel = 0; iel < n_elems_in_dom; iel++)
//     // {
//     //     int globEl = elemInDom[iel];
//     //      for (int iloc = 0; iloc < dim+1; iloc++)
//     //     {
//     //         int iglob = elem_v[globEl][iloc];
//     //         int optNode = optNodeFromGlobNode[iglob];
//     //         newDf0[optNode] += df0[optNode]*Volume_v[globEl]/(dim+1);
//     //         countW[optNode] += Volume_v[globEl]/(dim+1);
//     //     }
//     // }
//     // newDf0 /= countW;
//     // df0 = newDf0;
//     //--------------------------------------------------------------
//     // CONSTRAINTS DEFINITION
//     //--------------------------------------------------------------
//     // volume constraint
//     //--------------------
//     // int_\Omega { x } - V_r*V_0 <= 0 
//     prec integral = 0;
//     for (int iel = 0; iel < n_elems_in_dom; iel++)
//     {
//         int globEl = elemInDom[iel];
//         for (int iloc = 0; iloc < dim+1; iloc++)
//         {
//             int iglob = elem_v[globEl][iloc];
//             int optNode = optNodeFromGlobNode[iglob];
//             integral += x[optNode]/(dim+1)*Volume_v[globEl];
//         }
//     }
    
//     Vol = integral;
//     g[0]  = (integral - Vr * V0)/V0;
//     //--------------------------------------------------------------
//     // CONSTRAINTS DERIVATIVE
//     //--------------------------------------------------------------
//     dg.resetZeros();

//     for (int iel = 0; iel < n_elems_in_dom; iel++)
//     {
//         int globEl = elemInDom[iel];
//         for (int iloc = 0; iloc < dim+1; iloc++)
//         {
//             int iglob = elem_v[globEl][iloc];
//             int optNode = optNodeFromGlobNode[iglob];
//             dg[0][optNode] += Volume_v[globEl]/(dim+1)/V0;
//         }
//     }
// }

// void OPTIMIZER::update_val_and_derivative(VECTOR &x, prec &f0, VECTOR &df0, VECTOR &g, MATRIX& dg, prec &Vol)
// {
//     //pause();
//     obj_functional = 0.0;
//     func_in_box = 0.0;
//     f0 = 0.0;
    
//     func_val.resetZeros();
//     no_weighted_func_val.resetZeros();

//     d_obj_functional.resetZeros();
//     df0.resetZeros();

//     gamma = x; // copy the optimzation procedure value in gamma

//     topology_optimization_diffusive_filter();

//     eval_gamma_acc_and_derivative();

//     x = gamma; // copy the filtered and projected value of gamma into the optimization procedure (called gammaOpt in TopOpt.cpp)

//     int nElem_v = (*physics).nElem_v;
//     int nNode_v = (*physics).nNodes_v;
//     int n_node   = (*physics).nNodes;
//     int n_in_bd_elem = (*physics).n_inlet_bounds_elems;

//     MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
//     VECTOR Volume_v = (*physics).Volume_v;
//     MATRIX_INT inlet_bound_elem(n_in_bd_elem, dim, (*physics).inlet_bounds_elems.PP, (*physics).inlet_bounds_elems.P);
//     VECTOR area_inlet_bound_elem((*physics).area_inlet_bounds_elems.P, n_in_bd_elem);
    
//     VECTOR sol_times = (*physics).solution_times;
//     int n_times = sol_times.length;
//     VECTOR time_steps = (*physics).solution_deltaT;
//     prec tot_time = abs(sol_times.get_last() - sol_times[0]);
//     VECTOR time_weights(n_times);
//     if (time_integration_procedure == 0) // RECTANGLES RULE OVER TIME: only one suitable for stationary cases
//     {
//         time_weights[0] = 0;
//         for (int itime = 1; itime < n_times; itime++)
//         {
//             time_weights[itime] = time_steps[itime-1];
//         } 
//     }
//     else if (time_integration_procedure == 1) // TRAPEZIODAL RULE OVER TIME: not suitable for stationary problems, but it can be more reliable on time dependent ones
//     {
//         time_weights[0] = time_steps[0] / 2;
//         for (int itime = 1; itime < (n_times-1); itime++)
//         {
//             time_weights[itime] = (time_steps[itime-1] + time_steps[itime]) / 2;
//         } 
//         time_weights[n_times-1] = time_steps.get_last() / 2;
//     }
    
//     for (int itime = 0; itime < n_times; itime++)
//     {
//         prec curr_time_weigth = time_weights[itime];

//         VECTOR temp_NS_sol = (*physics).NS_solution.get_row(itime);
//         VECTOR temp_ADJ_sol = (*physics).NS_solution.get_row(itime);

//         MATRIX U(dim, nNode_v);
//         VECTOR P(n_node);
//         MATRIX Ua(dim, nNode_v);
//         VECTOR Pa(n_node);

//         decompose_solution(temp_NS_sol, U, P);
//         decompose_solution(temp_ADJ_sol, Ua, Pa);

//         temp_obj_functional = 0.0;
//         temp_func_val.resetZeros();
//         temp_d_obj_functional.resetZeros();

//         // get_functional_and_opt_derivative(elem_v, Volume_v, U, Ua, inlet_bound_elem, area_inlet_bound_elem, P);
//         eval_alpha_and_dAlpha();
//         eval_functional(elem_v, Volume_v, U, inlet_bound_elem, area_inlet_bound_elem, P);
//         eval_opt_derivative(elem_v, Volume_v, U, Ua);

//         obj_functional += temp_obj_functional * curr_time_weigth; 
//         func_in_box += temp_func_in_box * curr_time_weigth;
//         d_obj_functional += temp_d_obj_functional * curr_time_weigth;
//         func_val += (temp_func_val * curr_time_weigth);
//         no_weighted_func_val += (temp_no_weighted_func_val * curr_time_weigth);
//     }
//     obj_functional /= tot_time; 
//     func_in_box /= tot_time;
//     d_obj_functional /= tot_time;
//     func_val = func_val / tot_time;
    
//     if (first_it_flag == 1)
//     {
//         f0Init = obj_functional;
//         func_val_init = func_val;
//         func_in_box_init = func_in_box;
//         func_out_box_init = f0Init - func_in_box_init;
//         if (abs(f0Init) < 1e-10)
//         {
//             functional_normalization_factor = 1;
//         } 
//         else 
//         {
//             functional_normalization_factor = abs(f0Init);
//         }
//         (*physics).f0Init = f0Init;
//         (*physics).func_normalization_factor = functional_normalization_factor;
//         first_it_flag = 0;
//     }
    
//     obj_abs_val = obj_functional;
//     obj_functional /= functional_normalization_factor; 
//     func_in_box /= functional_normalization_factor;
//     func_out_box = obj_functional - func_in_box;
//     func_val = func_val / functional_normalization_factor;
//     d_obj_functional /= functional_normalization_factor;

//     f0 = obj_functional;
//     for (int ioptnode = 0; ioptnode < n_nodes_in_dom; ioptnode++)
//     {
//         int iglob = nodeInDom[ioptnode];
//         df0[ioptnode] = d_obj_functional[iglob];
//     }

//     //--------------------------------------------------------------
//     // CONSTRAINTS DEFINITION
//     //--------------------------------------------------------------
//     // volume constraint
//     //--------------------
//     // int_\Omega { gamma } - V_r*V_0 <= 0 
//     prec integral = 0;
//     for (int iel = 0; iel < n_elems_in_dom; iel++)
//     {
//         int globEl = elemInDom[iel];
//         for (int iloc = 0; iloc < dim+1; iloc++)
//         {
//             int iglob = elem_v[globEl][iloc];
//             int optNode = optNodeFromGlobNode[iglob];
//             integral += gamma_acc[optNode]/(dim+1)*Volume_v[globEl];
//         }
//     }
//     Vol = integral;
//     g[0]  = (integral - Vr * V0)/V0;
    
//     //--------------------------------------------------------------
//     // CONSTRAINTS DERIVATIVE
//     //--------------------------------------------------------------
//     dg.resetZeros();

//     for (int iel = 0; iel < n_elems_in_dom; iel++)
//     {
//         int globEl = elemInDom[iel];
//         for (int iloc = 0; iloc < dim+1; iloc++)
//         {
//             int iglob = elem_v[globEl][iloc];
//             int optNode = optNodeFromGlobNode[iglob];
//             dg[0][optNode] += dgamma_acc[optNode] * Volume_v[globEl]/(dim+1)/V0;
//         }
//     }
// }

void OPTIMIZER::update_val_and_derivative(VECTOR &x, prec &f0, VECTOR &df0, VECTOR &g, MATRIX& dg, prec &Vol)
{
   
    //--------------------------------------------------------------
    // EVALUATE FUNCTIONAL (and its derivatives) AND CONSTRAINTS
    //--------------------------------------------------------------
    update_val(x, f0, g, Vol);

    df0.resetZeros();
    //--------------------------------------------------------------
    // FUNCTIONAL DERIVATIVES
    // The derivatives are evaluated in the update_val method (even if it is not always necessary), 
    // so it remains just to copy their values also in the correct optimization enetiy: df0.
    // This procedure is redundand and not optimal, but allows to have a clear code in a complex part at a very small time cost. (that was the purpose at the time of the implementation)
    //--------------------------------------------------------------
    for (int ioptnode = 0; ioptnode < n_nodes_in_dom; ioptnode++)
    {
        int iglob = nodeInDom[ioptnode];
        df0[ioptnode] = d_obj_functional[iglob];
    }
    
    //--------------------------------------------------------------
    // CONSTRAINTS DERIVATIVES
    //--------------------------------------------------------------
    update_constraints_derivative(dg);
}

void OPTIMIZER::update_val(VECTOR &x, prec &f0, VECTOR &g, prec &Vol)
{
    obj_functional = 0.0;
    func_in_box = 0.0;
    f0 = 0.0;
    
    func_val.resetZeros();
    no_weighted_func_val.resetZeros();

    d_obj_functional.resetZeros();

    gamma = x; // copy the optimzation procedure value in gamma

    topology_optimization_diffusive_filter();

    eval_gamma_acc_and_derivative();
    // gamma_acc.printRowMatlab("acc");
    // pause();

    x = gamma_acc; // copy the filtered and projected value of gamma into the optimization procedure (called gammaOpt in TopOpt.cpp)

    MATRIX_INT elem_v((*physics).nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
    VECTOR Volume_v = (*physics).Volume_v;
    MATRIX_INT inlet_bound_elem((*physics).n_inlet_bounds_elems, dim, (*physics).inlet_bounds_elems.PP, (*physics).inlet_bounds_elems.P);
    VECTOR area_inlet_bound_elem((*physics).area_inlet_bounds_elems.P, (*physics).n_inlet_bounds_elems);

    switch ((*physics).isStationary)
    {
        case 1: // stationary solution
        {
            get_time_functional(0, elem_v, Volume_v, inlet_bound_elem, area_inlet_bound_elem);

            obj_functional = temp_obj_functional; 
            func_in_box = temp_func_in_box;
            d_obj_functional = temp_d_obj_functional;
            func_val = temp_func_val;
            no_weighted_func_val = temp_no_weighted_func_val;
            break;
        }
        default: // time dependent solution
        {
            VECTOR sol_times = (*physics).solution_times;
            int n_times = sol_times.length;
            VECTOR time_steps = (*physics).solution_deltaT;
            prec tot_time = abs(sol_times.get_last() - sol_times[0]);
            VECTOR time_weights(n_times);
            if (time_integration_procedure == 0) // RECTANGLES RULE OVER TIME: only one suitable to simulate the behavior of stationary cases
            {
                time_weights[0] = 0;
                for (int itime = 1; itime < n_times; itime++)
                {
                    time_weights[itime] = time_steps[itime-1];
                } 
            }
            else if (time_integration_procedure == 1) // TRAPEZIODAL RULE OVER TIME: not suitable for stationary problems, but it can be more reliable on time dependent ones
            {
                time_weights[0] = time_steps[0] / 2;
                for (int itime = 1; itime < (n_times-1); itime++)
                {
                    time_weights[itime] = (time_steps[itime-1] + time_steps[itime]) / 2;
                } 
                time_weights[n_times-1] = time_steps.get_last() / 2;
            }
            
            for (int itime = 0; itime < n_times; itime++)
            {
                prec curr_time_weigth = time_weights[itime];

                get_time_functional(itime, elem_v, Volume_v, inlet_bound_elem, area_inlet_bound_elem);

                // VECTOR temp_NS_sol = (*physics).NS_solution.get_row(itime);
                // VECTOR temp_ADJ_sol = (*physics).NS_solution.get_row(itime);

                // MATRIX U(dim, nNode_v);
                // VECTOR P(n_node);
                // MATRIX Ua(dim, nNode_v);
                // VECTOR Pa(n_node);

                // decompose_solution(temp_NS_sol, U, P);
                // decompose_solution(temp_ADJ_sol, Ua, Pa);

                // temp_obj_functional = 0.0;
                // temp_func_val.resetZeros();
                // temp_d_obj_functional.resetZeros();

                // get_functional_and_opt_derivative(elem_v, Volume_v, U, Ua, inlet_bound_elem, area_inlet_bound_elem, P);
                // eval_alpha_and_dAlpha();
                // eval_functional(elem_v, Volume_v, U, inlet_bound_elem, area_inlet_bound_elem, P);
                // eval_opt_derivative(elem_v, Volume_v, U, Ua);

                obj_functional += (temp_obj_functional * curr_time_weigth); 
                func_in_box += (temp_func_in_box * curr_time_weigth);
                d_obj_functional += (temp_d_obj_functional * curr_time_weigth);
                func_val += (temp_func_val * curr_time_weigth);
                no_weighted_func_val += (temp_no_weighted_func_val * curr_time_weigth);
            }
            obj_functional /= tot_time; 
            func_in_box /= tot_time;
            d_obj_functional /= tot_time;
            func_val = func_val / tot_time;
            break;
        }   
    }
    
    if (first_it_flag == 1)
    {
        f0Init = obj_functional;
        func_val_init = func_val;
        func_in_box_init = func_in_box;
        func_out_box_init = f0Init - func_in_box_init;
        if (abs(f0Init) < 1e-10)
        {
            functional_normalization_factor = 1;
        } 
        else 
        {
            functional_normalization_factor = abs(f0Init);
        }
        (*physics).f0Init = f0Init;
        (*physics).func_normalization_factor = functional_normalization_factor;
        first_it_flag = 0;
    }
    
    obj_abs_val = obj_functional;
    obj_functional /= functional_normalization_factor; 
    func_in_box /= functional_normalization_factor;
    func_out_box = obj_functional - func_in_box;
    func_val = func_val / functional_normalization_factor;
    d_obj_functional /= functional_normalization_factor;

    f0 = obj_functional;
    // for (int ioptnode = 0; ioptnode < n_nodes_in_dom; ioptnode++)
    // {
    //     int iglob = nodeInDom[ioptnode];
    //     df0[ioptnode] = d_obj_functional[iglob];
    // }

    //--------------------------------------------------------------
    // CONSTRAINTS EVALUATION
    //--------------------------------------------------------------
    update_constraints(g, elem_v, Volume_v);
    
    Vol = constraints.list[0].vol;
    
    // //--------------------------------------------------------------
    // // CONSTRAINTS DERIVATIVE
    // //--------------------------------------------------------------
    // dg.resetZeros();

    // for (int iel = 0; iel < n_elems_in_dom; iel++)
    // {
    //     int globEl = elemInDom[iel];
    //     for (int iloc = 0; iloc < dim+1; iloc++)
    //     {
    //         int iglob = elem_v[globEl][iloc];
    //         int optNode = optNodeFromGlobNode[iglob];
    //         dg[0][optNode] += dgamma_acc[optNode] * Volume_v[globEl]/(dim+1)/V0;
    //     }
    // }
}

void OPTIMIZER::update_constraints(VECTOR &g, MATRIX_INT &elem_v, VECTOR &Volume_v)
{
    for (int icons = 0; icons < constraints.n_constr; icons++)
    {    
        int type = constraints.list[icons].type;
        switch (type)
        {
            case 0:
            {
                update_volume_constraint(g, icons, constraints.list[icons], elem_v, Volume_v);
                break;
            }
            case 1:
            {
                update_subdomain_volume_constraint(g, icons, constraints.list[icons], elem_v, Volume_v);
                break;
            }
            case 2:
            {
                update_edge_size_constraint(g, icons, constraints.list[icons]);
                break;
            }
            case 3:
            {
                update_discretizing_constraint(g, icons, constraints.list[icons], elem_v, Volume_v);
                break;
            }
            default:
            {
                throw_line("ERROR: not handled constraint case\n");
                break;
            }  
        }
    }
}

void OPTIMIZER::update_constraints_derivative(MATRIX &dg)
{
    dg.resetZeros();
    for (int icons = 0; icons < constraints.n_constr; icons++)
    {    
        int type = constraints.list[icons].type;
        switch (type)
        {
            case 0:
            {
                update_volume_constraint_derivative(dg, icons, constraints.list[icons]);
                break;
            }
            case 1:
            {
                update_subdomain_volume_constraint_derivative(dg, icons, constraints.list[icons]);
                break;
            }
            case 2:
            {
                update_edge_size_constraint_derivative(dg, icons, constraints.list[icons]);
                break;
            }
            case 3:
            {
                update_discretizing_constraint_derivative(dg, icons, constraints.list[icons]);
                break;
            }
            default:
            {
                break;
            }  
        }
    }
}

void OPTIMIZER::update_volume_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr, MATRIX_INT &elem_v, VECTOR &Volume_v)
{
    // volume constraint
    //--------------------
    // int_\Omega { gamma } - V_r*V_0 <= 0 
    prec integral = 0;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            int optNode = optNodeFromGlobNode[iglob];
            integral += gamma_acc[optNode]/(dim+1)*Volume_v[globEl];
        }
    }
    constr.vol = integral;
    if (constr.vol_0 < 0.0)
    {
        constr.vol_0 = V0;
    }
    Vr = constr.Vr;
    g[iconstr]  = (integral - Vr * V0)/V0;
}

void OPTIMIZER::update_subdomain_volume_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr, MATRIX_INT &elem_v, VECTOR &Volume_v)
{
    // subdomain volume constraint
    //--------------------
    // int_\Omega_i { gamma } - V_r*V_0 <= 0 
    VECTOR_INT glob_elems_in_subdom = (*physics).elems_in_doms[constr.domain_id];
    prec integral = 0;
    for (int iel = 0; iel < glob_elems_in_subdom.length; iel++)
    {
        int globEl = glob_elems_in_subdom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            int optNode = optNodeFromGlobNode[iglob];
            integral += gamma_acc[optNode]/(dim+1)*Volume_v[globEl];
        }
    }
    constr.vol = integral;
    if (constr.vol_0 < 0.0)
    {
        constr.vol_0 = 0;
        for (int iel = 0; iel < glob_elems_in_subdom.length; iel++)
        {
            int globEl = glob_elems_in_subdom[iel];
            constr.vol_0 += Volume_v[globEl];
        }
    }
    g[iconstr]  = (integral - constr.Vr * constr.vol_0)/constr.vol_0;
}

void OPTIMIZER::update_discretizing_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr, MATRIX_INT &elem_v, VECTOR &Volume_v)
{
    // volume constraint
    //--------------------
    // int_\Omega { gamma } - V_r*V_0 <= 0 
    prec integral = 0;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            int optNode = optNodeFromGlobNode[iglob];
            prec temp_gamma = gamma_acc[optNode];
            integral += temp_gamma*(1-temp_gamma)/(dim+1)*Volume_v[globEl];
        }
    }

    constr.discretization_res = (integral / (*physics).V0) - constr.discretization_toll;
    g[iconstr]  = constr.discretization_res;
}

void OPTIMIZER::update_volume_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr)
{
    MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
    VECTOR Volume_v = (*physics).Volume_v;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            int optNode = optNodeFromGlobNode[iglob];
            dg[iconstr][optNode] += dgamma_acc[optNode] * Volume_v[globEl]/(dim+1)/V0;
        }
    }
}

void OPTIMIZER::update_subdomain_volume_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr)
{
    MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
    VECTOR Volume_v = (*physics).Volume_v;
    VECTOR_INT glob_elems_in_subdom = (*physics).elems_in_doms[constr.domain_id];
    prec vol_0 = constr.vol_0;
    for (int iel = 0; iel < glob_elems_in_subdom.length; iel++)
    {
        int globEl = glob_elems_in_subdom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            int optNode = optNodeFromGlobNode[iglob];
            dg[iconstr][optNode] += dgamma_acc[optNode] * Volume_v[globEl]/(dim+1)/vol_0;
        }
    }
}

void OPTIMIZER::update_discretizing_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr)
{
    MATRIX_INT elem_v(nElem_v, dim+1, (*physics).elem_v.PP, (*physics).elem_v.P);
    VECTOR Volume_v = (*physics).Volume_v;
    prec vol_0 = (*physics).V0;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            int optNode = optNodeFromGlobNode[iglob];
            prec temp_coef = 1 - 2*gamma_acc[optNode];
            dg[iconstr][optNode] += temp_coef * dgamma_acc[optNode] * Volume_v[globEl]/(dim+1)/vol_0;
        }
    }
}

void OPTIMIZER::update_edge_size_constraint(VECTOR &g, int iconstr, CONSTRAINT &constr)
{
    int bound_id = constr.bound_id;
    MATRIX_INT bound_elems = (*physics).bounds_elems_v[bound_id];
    int n_el_in_bound = bound_elems.nRow;
    VECTOR elem_gamma_avg(n_el_in_bound);
    for (int iel = 0; iel < n_el_in_bound; iel++)
    {
        prec temp_gamma_avg = 0;
        for (int inode = 0; inode < dim; inode++)
        {
            int opt_node = optNodeFromGlobNode[bound_elems[iel][inode]];
            temp_gamma_avg += gamma_acc[opt_node];
        }
        temp_gamma_avg /= dim;
        elem_gamma_avg[iel] = temp_gamma_avg;
    }
    prec integral = VECTOR::dot(elem_gamma_avg.P, (*physics).bounds_elems_surface_v[bound_id].P, n_el_in_bound);
    constr.surf = integral;
    if (constr.surf_0 < 0.0)
    {
        constr.surf_0 = (*physics).bounds_elems_surface_v[bound_id].sum();
    }
    prec S0 = constr.surf_0;
    prec Sr = constr.Sr;
    g[iconstr] = (integral - Sr * S0)/S0;
}

void OPTIMIZER::update_edge_size_constraint_derivative(MATRIX &dg, int iconstr, CONSTRAINT &constr)
{
    int bound_id = constr.bound_id;
    prec S0 = constr.surf_0;
    MATRIX_INT bound_elems = (*physics).bounds_elems_v[bound_id];
    int n_el_in_bound = bound_elems.nRow;
    VECTOR elem_gamma_avg(n_el_in_bound);
    for (int iel = 0; iel < n_el_in_bound; iel++)
    {
        for (int inode = 0; inode < dim; inode++)
        {
            int opt_node = optNodeFromGlobNode[bound_elems[iel][inode]];
            dg[iconstr][opt_node] += dgamma_acc[opt_node] * (*physics).bounds_elems_surface_v[bound_id][iel] / dim / S0;
        }
    }
}

void OPTIMIZER::decompose_solution(VECTOR &sol, MATRIX &U_sol, VECTOR &P_sol)
{
    int nNnodes_v = (*physics).nNodes_v;
    int nNodes    = (*physics).nNodes;
    for (int icomp = 0; icomp < dim; icomp++)
    {
        int start_v_comp_id = icomp * nNodes_v;
        for (int inode = 0; inode < nNodes_v; inode++)
        {
            U_sol[icomp][inode] = sol[start_v_comp_id + inode];
        }
    }
    int start_p_sol_id = dim*nNnodes_v;
    for (int inode = 0; inode < nNodes; inode++)
    {
        P_sol[inode] = sol[start_p_sol_id + inode];
    }
}

//-------------------------------------------------------
// GOC
//-------------------------------------------------------
void OPTIMIZER::solveGOC(VECTOR &x, prec &Vol, prec &f0)
{
    std::cout << "\n-----| SOLVE GOC |-----\n";
    static int nCons = 1;
    static VECTOR lambda = VECTOR::zeros(nCons) + 1;
    static VECTOR gold = VECTOR::zeros(nCons);
    static VECTOR df0(n_nodes_in_dom);
    static VECTOR g(nCons);
    static MATRIX dg(nCons, n_nodes_in_dom);

    update_val_and_derivative(x, f0, df0, g, dg, Vol);

    prec eps = 0.05; prec maxmove = 0.2;

    VECTOR deltag = g - gold; gold = g;
    VECTOR xNew(n_nodes_in_dom);
    //------------------------
    for (int icons = 0; icons < nCons; icons++)
    {
        prec p0;
        if ((g[icons] > 0 && deltag[icons] > 0)  || (g[icons] < 0 && deltag[icons] < 0)) p0 = 1.0;
        else if ( (g[icons] > 0 && deltag[icons] > -eps) || (g[icons] < 0 && deltag[icons] < eps)) p0 = 0.5;
        else p0 = 0;
        lambda[icons] *= (1+ p0*(g[icons] + deltag[icons]));
    }
    //-------------------------------------------------------
    for (int iglob = 0; iglob < n_nodes_in_dom; iglob++)
    {
        prec num = 0;
        prec den = 0;
        if (df0[iglob] < 0) num = df0[iglob];
        else den = df0[iglob];

        for (int icons = 0; icons < nCons; icons++)
        {
            if (dg[icons][iglob] < 0) num += lambda[icons] * dg[icons][iglob];
            else den += lambda[icons] * dg[icons][iglob];
        }

        prec De = - num/den;
        xNew[iglob] = max(1e-6,max(x[iglob] - maxmove, min(1.0,min(x[iglob] + maxmove, x[iglob]*sqrt(De)))));
    }
    x = xNew;

    update_val_and_derivative(x, f0, df0, g, dg, Vol);
}

//-------------------------------
// !!!!!!!!!!!!! MMA !!!!!!!!!!!!
//-------------------------------
void OPTIMIZER::solveMMA(VECTOR &x, prec &Vol, prec &f0)
{
    std::cout << "\n-----| SOLVE MMA |-----\n";
    int nCons = constraints.n_constr;
    VECTOR df0(n_nodes_in_dom);
    VECTOR g;
    g.setZeros(nCons);
    MATRIX dg;
    dg.zeros(nCons, n_nodes_in_dom);
    //---
    int m = nCons;
    int n = n_nodes_in_dom;
    // prec epsimin = 0.0000001;
    VECTOR xval;
    xval    = x;
    VECTOR xold1;
    xold1   = xval;
    VECTOR xold2;
    xold2   = xval;
    VECTOR xmin    = VECTOR::zeros(n_nodes_in_dom);
    VECTOR xmax    = VECTOR::zeros(n_nodes_in_dom) + 1;
    VECTOR low;
    low     = xmin;
    VECTOR upp;
    upp     = xmax;
    VECTOR c(nCons); //c[0] = 1000;
    VECTOR d(nCons); //d[0] = 1;
    int a0       = 1;
    VECTOR a(nCons); //a[0] = 0;
    for (int i = 0; i < nCons; i++)
    {
        c[i] = 1000;
        d[i] = 1;
        a[i] = 0;
    }
    int outeriter = 0;
    int maxoutit  = 1;
    int kkttol    = 1e-12;
    //---
    if (outeriter < 0.5)
    {
        update_val_and_derivative(x, f0, df0, g, dg, Vol);
    }
    prec kktnorm = kkttol+10;
    int outit = 0;
    //---
    VECTOR xmma(n_nodes_in_dom);  VECTOR ymma(nCons); prec zmma;
    VECTOR lam(nCons); VECTOR eta; prec zet; VECTOR mu(nCons); VECTOR s(nCons);
    VECTOR xsi(n_nodes_in_dom);
    while (kktnorm > kkttol && outit < maxoutit)
    {
        outit++;
        outeriter++;
        // The MMA subproblem is solved at the point xval:
        mmasub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp,
               m,n,outeriter,xval,xmin,xmax,xold1,xold2, f0,df0,g,dg,a0,a,c,d);
        xold2 = xold1;
        xold1 = xval;
        xval = xmma;
        update_val_and_derivative(xval, f0, df0, g, dg, Vol);
        kktnorm = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0,g,dg,a0,a,c,d);
    }
    x = xval;
    // free memory
}

//----------------------------------------
void OPTIMIZER::mmasub(VECTOR &xmma, VECTOR &ymma, prec &zmma, VECTOR &lam, VECTOR &xsi, VECTOR &eta, VECTOR &mu, prec &zet, VECTOR &s, VECTOR &low, VECTOR &upp,
                int m,int n,int iter, VECTOR xval,VECTOR xmin,VECTOR xmax,VECTOR xold1,VECTOR xold2, prec f0, 
                VECTOR df0, VECTOR g, MATRIX dg, prec a0, VECTOR a, VECTOR c,VECTOR d)
{
    prec epsimin = 1e-7;
    prec raa0 = 0.00001;
    prec move = 0.5;
    prec albefa = 0.1;
    prec asyinit = 0.5;
    prec asyincr = 1.2;
    prec asydecr = 0.7;

    VECTOR eeen = VECTOR::zeros(n) + 1;
    VECTOR eeem = VECTOR::zeros(m) + 1;
    VECTOR zeron = VECTOR::zeros(n);

    VECTOR zzz(n);
    // Calculation of the asymptotes low and upp :
    if (iter < 2.5)
    {
        VECTOR tempVec = xmax-xmin; tempVec *= asyinit;
        low = xval - tempVec;
        upp = xval + tempVec;
    }
    else
    {
        zzz = (xval-xold1)*(xold1-xold2);
        VECTOR factor = eeen;

        for (int i = 0; i < n; i++)
        {
            if (zzz[i] > 0) factor[i] = asyincr;
            else if (zzz[i] < 0) factor[i] = asydecr;
        }
        low = xval - factor*(xold1 - low);
        upp = xval + factor*(upp - xold1);
        VECTOR lowmin = xval - (xmax-xmin)*10;
        VECTOR lowmax = xval - (xmax-xmin)*0.01;
        VECTOR uppmin = xval + (xmax-xmin)*0.01;
        VECTOR uppmax = xval + (xmax-xmin)*10;
        low = max(low,lowmin);
        low = min(low,lowmax);
        upp = min(upp,uppmax);
        upp = max(upp,uppmin);
    }
    
    // Calculation of the bounds alfa and beta :

    VECTOR zzz1;
    zzz1 = low;
    VECTOR temp = (xval-low); temp *= albefa;
    zzz1 += temp;
    VECTOR zzz2 = xval - (xmax-xmin)*move;
    zzz  = max(zzz1,zzz2);
    VECTOR alfa = max(zzz,xmin);
    zzz1 = upp - (upp-xval)*albefa;
    zzz2 = xval + (xmax-xmin)*move;
    zzz  = min(zzz1,zzz2);
    VECTOR beta = min(zzz,xmax);
    
    // Calculations of p0, q0, P, Q and b.
    VECTOR xmami = xmax-xmin;
    VECTOR xmamieps = eeen*0.00001;
    xmami = max(xmami,xmamieps);
    VECTOR xmamiinv = eeen.pointdiv(xmami);
    VECTOR ux1 = upp-xval;
    VECTOR ux2 = ux1*ux1;
    VECTOR xl1 = xval-low;
    VECTOR xl2 = xl1*xl1;
    VECTOR uxinv = eeen.pointdiv(ux1);
    VECTOR xlinv = eeen.pointdiv(xl1);
    
    //
    VECTOR p0;
    p0 = zeron;
    VECTOR q0;
    q0 = zeron;
    p0 = max(df0,0);
    q0 = max(df0*(-1),0);
    VECTOR pq0 = (p0 + q0)*0.001 + xmamiinv*raa0;
    p0 = p0 + pq0;
    q0 = q0 + pq0;
    p0 = p0*ux2;
    q0 = q0*xl2;
    //
    int nCons = dg.nRow;
    std::vector<VECTOR> P(nCons);
    for (int i = 0; i < nCons; i++)
    {
        P[i] = max(dg[i], 0, n);
    }
    std::vector<VECTOR> Q(nCons);
    for (int i = 0; i < nCons; i++)
    {

        Q[i].initialize(n);
        Q[i] = min(dg[i], 0, n);
        Q[i] *= -1;
    }
    VECTOR b = VECTOR::zeros(nCons);
    for (int i = 0; i < nCons; i++)
    {
        VECTOR PQ = (P[i] + Q[i])*0.001 + xmamiinv*raa0;
        
        P[i] += PQ; Q[i] += PQ;
        P[i] *= ux2;
        Q[i] *= xl2;

        b[i] += P[i].dot(uxinv) + Q[i].dot(xlinv);
    }
    b -= g;

    // Solving the subproblem by a primal-dual Newton method
    subsolv(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,
            m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);

}


void OPTIMIZER::subsolv(VECTOR &xmma, VECTOR &ymma, prec &zmma, VECTOR &lamma, VECTOR &xsimma, VECTOR &etamma, VECTOR &mumma, prec &zetmma, VECTOR &smma,
                int m,int n, prec epsimin, VECTOR low, VECTOR upp, VECTOR alfa,VECTOR beta, VECTOR p0,VECTOR q0,std::vector<VECTOR> P,std::vector<VECTOR> Q, 
                prec a0,VECTOR a, VECTOR b,VECTOR c,VECTOR d)
{
    VECTOR een = VECTOR::zeros(n) + 1;
    VECTOR eem = VECTOR::zeros(m) + 1;
    prec epsi = 1;
    VECTOR epsvecn = een*epsi;
    VECTOR epsvecm = eem*epsi;
    VECTOR x = (alfa+beta)*0.5;
    VECTOR y;
    y = eem;
    prec z = 1;
    VECTOR lam;
    lam = eem;
    VECTOR xsi = een.pointdiv(x-alfa);
    xsi = max(xsi,een);
    VECTOR eta = een.pointdiv(beta-x);
    eta = max(eta,een);
    VECTOR mu  = max(eem,c*0.5);
    prec zet = 1;
    VECTOR s;
    s = eem;
    int itera = 0;

    while (epsi > epsimin)
    {
        epsvecn = een*epsi;
        epsvecm = eem*epsi;
        VECTOR ux1 = upp-x;
        VECTOR xl1 = x-low;
        VECTOR ux2 = ux1*ux1;
        VECTOR xl2 = xl1*xl1;
        VECTOR uxinv1 = een.pointdiv(ux1);
        VECTOR xlinv1 = een.pointdiv(xl1);
        VECTOR plam;
        plam = p0;
        VECTOR qlam;
        qlam = q0;

        VECTOR gvec(m);
        for (int j = 0; j < m; j++)
        {
            VECTOR tempP;
            tempP = P[j];
            VECTOR tempQ;
            tempQ = Q[j];
            plam += tempP*lam[j]; 
            qlam += tempQ*lam[j];

            gvec[j] = tempP.dot(uxinv1) + tempQ.dot(xlinv1);
        }
        VECTOR dpsidx = plam.pointdiv(ux2) - qlam.pointdiv(xl2);
        VECTOR rex = dpsidx - xsi + eta;
        VECTOR rey = c + d*y - mu - lam;
        prec rez = a0 - zet - a.dot(lam);
        VECTOR relam = gvec - a*z - y + s - b;
        VECTOR rexsi = xsi*(x-alfa) - epsvecn;
        VECTOR reeta = eta*(beta-x) - epsvecn;
        VECTOR remu = mu*y - epsvecm;
        prec rezet = zet*z - epsi;
        VECTOR res = lam*s - epsvecm;
        // copy into residu
        VECTOR residu(3*n+4*m+2);
        prec* tempP = &(residu[0]);
        rex.copyTo(tempP); tempP += n; 
        rey.copyTo(tempP); tempP += m; 
        *tempP = rez; tempP++;
        relam.copyTo(tempP); tempP += m; 
        rexsi.copyTo(tempP); tempP += n; 
        reeta.copyTo(tempP); tempP += n; 
        remu.copyTo(tempP); tempP += m; 
        *tempP = rezet; tempP++;
        res.copyTo(tempP); tempP += m;
        //
        prec residunorm = residu.norm();
        prec residumax = (residu.vecAbs()).max();
        int ittt = 0;
        while (residumax > 0.9*epsi && ittt < 200)
        {
            ittt++;
            itera++;
            ux1 = upp-x;
            xl1 = x-low;
            ux2 = ux1*ux1;
            xl2 = xl1*xl1;
            VECTOR ux3 = ux1*ux2;
            VECTOR xl3 = xl1*xl2;
            uxinv1 = een.pointdiv(ux1);
            xlinv1 = een.pointdiv(xl1);
            VECTOR uxinv2 = een.pointdiv(ux2);
            VECTOR xlinv2 = een.pointdiv(xl2);
            plam = p0;
            qlam = q0;
            std::vector<VECTOR> GG(m);
            for (int j = 0; j < m; j++)
            {
                VECTOR tempP;
                tempP = P[j];
                VECTOR tempQ;
                tempQ = Q[j];
                GG[j].initialize(n);
                plam += tempP*lam[j]; 
                qlam += tempQ*lam[j];

                gvec[j] = tempP.dot(uxinv1) + tempQ.dot(xlinv1);

                GG[j] = tempP*(uxinv2) - tempQ*(xlinv2);
            }
            dpsidx = plam.pointdiv(ux2) - qlam.pointdiv(xl2);
            VECTOR delx = dpsidx - epsvecn.pointdiv(x-alfa) + epsvecn.pointdiv(beta-x);
            VECTOR dely = c + d*y - lam - epsvecm.pointdiv(y);
            prec delz = a0 - a.dot(lam) - epsi/z;
            VECTOR dellam = gvec - a*z - y - b + epsvecm.pointdiv(lam);
            VECTOR diagx = plam.pointdiv(ux3) + qlam.pointdiv(xl3);
            diagx = diagx*2 + xsi.pointdiv(x-alfa) + eta.pointdiv(beta-x);
            VECTOR diagxinv = een.pointdiv(diagx);
            VECTOR diagy = d + mu.pointdiv(y);
            VECTOR diagyinv = eem.pointdiv(diagy);
            VECTOR diaglam = s.pointdiv(lam);
            VECTOR diaglamyi = diaglam+diagyinv;
            VECTOR dlam(m);
            VECTOR dx(n);
            prec dz;
            if (m < n)
            {
                VECTOR blam = dellam + dely.pointdiv(diagy);
                for (int j = 0; j < m; j++)
                {
                    blam[j] -= GG[j].dot(delx.pointdiv(diagx));
                }  
                // blam.print("../../advanced fluid/codes/01-06-22/12-05-2022/blam.txt");
                VECTOR bb;
                bb = blam;
                
                bb.append(delz);
                MATRIX AA(m+1, m+1);

                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        AA[i][j] = (GG[i]*diagxinv).dot(GG[j]);
                        if (i == j) AA[i][j] += diaglamyi[i];
                    }  
                    AA[i][m] = a[i];
                }  
                for (int j = 0; j < m; j++) AA[m][j] = a[j];
                AA[m][m] = -zet/z;

                MATRIX QQ(m+1, m+1); UPPER_TRIANGULAR RR(m+1);
                QR_DECOMPOSITION::reduced(AA, QQ, RR);
                MATRIX newQQ = QQ.transpose();
                QQ = newQQ;
                bb = QQ*bb;
                VECTOR solut(RR.N);
                RR.solveLS(bb, solut);

                VECTOR tempVec = VECTOR::zeros(n);
                for (int i = 0; i < m; i++)
                {
                    dlam[i] = solut[i];
                    tempVec += GG[i]*dlam[i];
                }
                dz = solut[m];
                dx = delx.pointdiv(diagx)*(-1) - tempVec.pointdiv(diagx);
            }
            else
            {
                throw_line("ERROR IN MMA, more constraints than nodes, not implemented mode.\n");
            }
            //
            VECTOR dy = dely.pointdiv(diagy)*(-1) + dlam.pointdiv(diagy);
            VECTOR dxsi = xsi*(-1) + epsvecn.pointdiv(x-alfa) - (xsi*dx).pointdiv(x-alfa);
            VECTOR deta = eta*(-1) + epsvecn.pointdiv(beta-x) + (eta*dx).pointdiv(beta-x);
            VECTOR dmu  = mu*(-1) + epsvecm.pointdiv(y) - (mu*dy).pointdiv(y);
            prec dzet = -zet + epsi/z - zet*dz/z;
            VECTOR ds   = s*(-1)+ epsvecm.pointdiv(lam) - (s*dlam).pointdiv(lam);
            // copy for xx
            VECTOR xx(4*m+2+2*n);
            tempP = &(xx[0]);
            y.copyTo(tempP); tempP += m; 
            *tempP = z; tempP++; 
            lam.copyTo(tempP); tempP += m; 
            xsi.copyTo(tempP); tempP += n; 
            eta.copyTo(tempP); tempP += n; 
            mu.copyTo(tempP); tempP += m; 
            *tempP = zet; tempP++;
            s.copyTo(tempP); tempP += m; 
            // copy dor dxx
            VECTOR dxx(4*m+2+2*n);
            tempP = &(dxx[0]);
            dy.copyTo(tempP); tempP += m; 
            *tempP = dz; tempP++; 
            dlam.copyTo(tempP); tempP += m; 
            dxsi.copyTo(tempP); tempP += n; 
            deta.copyTo(tempP); tempP += n; 
            dmu.copyTo(tempP); tempP += m; 
            *tempP = dzet; tempP++;
            ds.copyTo(tempP); tempP += m; 
            //
            VECTOR stepxx = dxx.pointdiv(xx)*(-1.01);
            prec stmxx  = stepxx.max();
            VECTOR stepalfa = dx.pointdiv(x-alfa)*(-1.01);
            prec stmalfa = stepalfa.max();
            VECTOR stepbeta = dx.pointdiv(beta-x)*(1.01);
            prec stmbeta = stepbeta.max();
            prec stmalbe  = max(stmalfa,stmbeta);
            prec stmalbexx = max(stmalbe,stmxx);
            prec stminv = max(stmalbexx,1);
            prec steg = 1/stminv;
            //
            VECTOR xold;
            xold   =   x;
            VECTOR yold;
            yold   =   y;
            prec zold   =   z;
            VECTOR lamold;
            lamold =  lam;
            VECTOR xsiold;
            xsiold =  xsi;
            VECTOR etaold;
            etaold =  eta;
            VECTOR muold;
            muold  =  mu;
            prec zetold =  zet;
            VECTOR sold;
            sold   =   s;
            //
            int itto = 0;
            prec resinew = 2*residunorm;
            while (resinew > residunorm && itto < 50)
            {
                itto++;
                x   =   xold + dx*steg;
                y   =   yold + dy*steg;
                z   =   zold + dz*steg;
                lam = lamold + dlam*steg;
                xsi = xsiold + dxsi*steg;
                eta = etaold + deta*steg;
                mu  = muold  + dmu*steg;
                zet = zetold + dzet*steg;
                s   =   sold + ds*steg;
                ux1 = upp-x;
                xl1 = x-low;
                ux2 = ux1*ux1;
                xl2 = xl1*xl1;
                uxinv1 = een.pointdiv(ux1);
                xlinv1 = een.pointdiv(xl1);

                plam = p0;
                qlam = q0;
                for (int j = 0; j < m; j++)
                {
                    VECTOR tempP;
                    tempP = P[j];
                    VECTOR tempQ;
                    tempQ = Q[j];
                    plam += tempP*lam[j]; 
                    qlam += tempQ*lam[j];

                    gvec[j] = tempP.dot(uxinv1) + tempQ.dot(xlinv1);
                }

                dpsidx = plam.pointdiv(ux2) - qlam.pointdiv(xl2) ;
                rex = dpsidx - xsi + eta;
                rey = c + d*y - mu - lam;
                rez = a0 - zet - a.dot(lam);
                relam = gvec - a*z - y + s - b;
                rexsi = xsi*(x-alfa) - epsvecn;
                reeta = eta*(beta-x) - epsvecn;
                remu = mu*y - epsvecm;
                rezet = zet*z - epsi;
                res = lam*s - epsvecm;
                 // copy into residu
                prec* tempP = &(residu[0]);
                rex.copyTo(tempP); tempP += n; 
                rey.copyTo(tempP); tempP += m; 
                *tempP = rez; tempP++;
                relam.copyTo(tempP); tempP += m; 
                rexsi.copyTo(tempP); tempP += n; 
                reeta.copyTo(tempP); tempP += n; 
                remu.copyTo(tempP); tempP += m; 
                *tempP = rezet; tempP++;
                res.copyTo(tempP); tempP += m;
                resinew = residu.norm();
                steg /= 2;
            }
            residunorm=resinew;
            residumax = (residu.vecAbs()).max();
            steg *= 2;
            
        }
        if (ittt > 198) printf("MMA WARNING: ittt > 198\n");
        epsi *= 0.1;
    }
    xmma   =   x;
    ymma   =   y;
    zmma   =   z;
    lamma =  lam;
    xsimma =  xsi;
    etamma =  eta;
    mumma  =  mu;
    zetmma =  zet;
    smma   =   s;
    // free memory 
}

prec OPTIMIZER::kktcheck(int m, int n, VECTOR x, VECTOR y, prec z, VECTOR lam, VECTOR xsi, VECTOR eta, VECTOR mu, prec zet,VECTOR s,
                 VECTOR xmin, VECTOR xmax, VECTOR df0, VECTOR g, MATRIX dg, prec a0, VECTOR a, VECTOR c, VECTOR d)
{
    VECTOR rex   = df0 - xsi + eta;
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < m; j++) rex[i] += dg[j][i]*lam[j];
    }
    VECTOR rey   = c + d*y - mu - lam;
    prec rez   = a0 - zet - a.dot(lam);
    VECTOR relam = g - a*z - y + s;
    VECTOR rexsi = xsi*(x-xmin);
    VECTOR reeta = eta*(xmax-x);
    VECTOR remu  = mu*y;
    prec rezet = zet*z;
    VECTOR res   = lam*s;
    //
    VECTOR residu(3*n+4*m+2);
    prec* tempP = &(residu[0]);
    rex.copyTo(tempP); tempP += n; 
    rey.copyTo(tempP); tempP += m; 
    *tempP = rez; tempP++;
    relam.copyTo(tempP); tempP += m; 
    rexsi.copyTo(tempP); tempP += n; 
    reeta.copyTo(tempP); tempP += n; 
    remu.copyTo(tempP); tempP += m; 
    *tempP = rezet; tempP++;
    res.copyTo(tempP); tempP += m;
    //
    prec residunorm = residu.norm();
    // prec residumax = residu.vecAbs().max();

    // free memory
    return residunorm;
}


//-------------------------------
// !!!!!!!!!!!!! GCMMA !!!!!!!!!!!!
//-------------------------------
void OPTIMIZER::solveGCMMA(VECTOR &x, prec &Vol, prec &f0)
{
    static int nCons = 1;
    static VECTOR df0(n_nodes_in_dom);
    static VECTOR g(nCons);
    static MATRIX dg(nCons, n_nodes_in_dom);
    //---
    int m = 1; // Number of Constraints
    int n = n_nodes_in_dom;
    prec epsimin = 0.0000001;
    VECTOR xval;
    xval    = x;
    VECTOR xold1;
    xold1   = xval;
    VECTOR xold2;
    xold2   = xval;
    VECTOR xmin    = VECTOR::zeros(n_nodes_in_dom);
    VECTOR xmax    = VECTOR::zeros(n_nodes_in_dom) + 1;
    VECTOR low;
    low     = xmin;
    VECTOR upp;
    upp     = xmax;
    VECTOR c(nCons); c[0] = 1000;
    VECTOR d(nCons); d[0] = 1;
    int a0       = 1;
    VECTOR a(nCons); a[0] = 0;
    //---------------
    prec raa0    = 0.01;
    VECTOR raa(1);    raa[0] = 0.01;
    prec raa0eps = 0.000001;
    VECTOR raaeps(1); raaeps[0] = 0.000001;
    //--------------
    int outeriter = 0;
    int maxoutit  = 1;
    int kkttol    = 1e-12;
    int innerit = 0;
    //---
    if (outeriter < 0.5) //always done since outeriter = 0
    {
        innerit = 0;
        update_val_and_derivative(x, f0, df0, g, dg, Vol);
    }
    prec kktnorm = kkttol+10;
    int outit = 0;
    //---
    VECTOR xmma(n_nodes_in_dom);  VECTOR ymma(nCons); prec zmma = 0;
    VECTOR lam(nCons); VECTOR eta; prec zet = 0; VECTOR mu(nCons); VECTOR s(nCons);
    VECTOR xsi(n_nodes_in_dom);
    //----------------
    // GOC SPECIAL PARAMETERS
    //---------------
    prec f0app;
    VECTOR fapp(m);
    prec f0valnew; // functional value
    VECTOR fvalnew(m); // constraints
    prec conserv;
    while (kktnorm > kkttol && outit < maxoutit)
    {
        outit++;
        outeriter++;
        //The parameters low, upp, raa0 and raa are calculated:
        asynmp(low,upp,raa0,raa,
              outeriter,n,xval,xold1,xold2,xmin,xmax, raa0eps,raaeps,df0,dg);
        // The MMA subproblem is solved at the point xval:
        gcmmasub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp,
                 m,n,outeriter,epsimin,xval,xmin,xmax, low, upp, raa0, raa, f0,df0,g,dg,a0,a,c,d);
        // The user should now calculate function values (no gradients)
        // of the objective- and constraint functions at the point xmma
        // ( = the optimal solution of the subproblem).
        // The results should be put in f0valnew and fvalnew.
        // updateJustVal(xmma, f0valnew, fvalnew);
        update_val(xmma, f0valnew, fvalnew, Vol);
    
        // It is checked if the approximations are conservative:
        concheck(conserv,
                 m,epsimin,f0app,f0valnew,fapp,fvalnew);
        // While the approximations are non-conservative (conserv=0),
        // repeated inner iterations are made:
        innerit = 0;

        if (conserv == 0)
        {
            while (conserv == 0 && innerit <= 15)
            {
                innerit++;
                // New values on the parameters raa0 and raa are calculated:
                raaupdate(raa0,raa,
                          xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew,f0app,fapp,raa0eps,raaeps,epsimin);

                // The GCMMA subproblem is solved with these new raa0 and raa:
                gcmmasub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp,
                         m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, raa0,raa,f0,df0,g,dg,a0,a,c,d);
                // The user should now calculate function values (no gradients)
                // of the objective- and constraint functions at the point xmma
                // ( = the optimal solution of the subproblem).
                // The results should be put in f0valnew and fvalnew.
                // updateJustVal(xmma, f0valnew, fvalnew);
                update_val(xmma, f0valnew, fvalnew, Vol);
                // It is checked if the approximations have become conservative:
                concheck(conserv,
                        m,epsimin,f0app,f0valnew,fapp,fvalnew);
            }
        }
        // No more inner iterations. Some vectors are updated:
        xold2 = xold1;
        xold1 = xval;
        xval  = xmma;
        // The user should now calculate function values and gradients
        // of the objective- and constraint functions at xval.
        // The results should be put in f0val, df0dx, fval and dfdx:
        update_val_and_derivative(xval, f0, df0, g, dg, Vol);

        // The residual vector of the KKT conditions is calculated:

        kktnorm = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, xmin,xmax,df0,g,dg,a0,a,c,d);
    }
    x = xval;
    // free memory
}


void OPTIMIZER::asynmp(VECTOR &low, VECTOR &upp, prec &raa0, VECTOR& raa,
              int outeriter,int n,VECTOR &xval,VECTOR &xold1,VECTOR &xold2,VECTOR &xmin,VECTOR &xmax, prec raa0eps,VECTOR &raaeps,VECTOR &df0,MATRIX &dg)
{
    //
    VECTOR eeen(n); eeen.reset(1.0);
    VECTOR xmami = xmax - xmin;
    VECTOR xmamieps = eeen *  0.00001;
    xmami = max(xmami,xmamieps);

    VECTOR temp1 = (df0.vecAbs());
    raa0 = temp1.dot(xmami);
    raa0 = max(raa0eps,(0.1/n)*raa0);   // WARNING !!!!!!!!
    raa  = (dg.absMat())*xmami;            // WARNING !!!!!!!!

    VECTOR temp2 = raa*(0.1/n);
    raa  = max(raaeps, temp2);
    if (outeriter < 2.5)
    {
        low = xval - xmami*0.5;
        upp = xval + xmami*0.5;
    }
    else
    {
        VECTOR xxx = (xval-xold1)*(xold1-xold2);
        VECTOR factor = eeen;

        for (int i = 0; i < n; i++)
        {
            if (xxx[i] == 0) continue;
            else if (xxx[i] > 0) factor[i] = 1.2;
            else factor[i] = 0.7;
            
        }
        low = xval - factor*(xold1 - low);
        upp = xval + factor*(upp - xold1);
        VECTOR lowmin = xval - xmami*10;
        VECTOR lowmax = xval - xmami*0.01;
        VECTOR uppmin = xval + xmami*0.01;
        VECTOR uppmax = xval + xmami*10;
        low = max(low,lowmin);
        low = min(low,lowmax);
        upp = min(upp,uppmax);
        upp = max(upp,uppmin);
    }  
}

void OPTIMIZER::gcmmasub(VECTOR &xmma,VECTOR &ymma,prec zmma,VECTOR &lam, VECTOR &xsi,VECTOR &eta,VECTOR &mu,prec zet,VECTOR &s, prec f0app, VECTOR &fapp,
                 int m,int n,int outeriter,prec epsimin,VECTOR &xval,VECTOR &xmin,VECTOR &xmax, VECTOR &low, VECTOR &upp, prec raa0, VECTOR &raa, prec f0,VECTOR &df0,
                 VECTOR &g,MATRIX &dg,prec a0,VECTOR &a,VECTOR & c, VECTOR &d)
{
    //
    VECTOR eeen(n); eeen.reset(1.0);
    VECTOR zeron(n); zeron.resetZeros();
    //
    // Calculations of the bounds alfa and beta.
    prec albefa = 0.1;
    prec move = 0.5;
    //
    VECTOR zzz1 = low + (xval-low)*albefa;
    VECTOR zzz2 = xval - (xmax-xmin)*move;
    VECTOR zzz  = max(zzz1,zzz2);
    VECTOR alfa = max(zzz,xmin);
    zzz1 = upp - (upp-xval)*albefa;
    zzz2 = xval + (xmax-xmin)*move;
    zzz  = min(zzz1,zzz2);
    VECTOR beta = min(zzz,xmax);
    //
    // Calculations of p0, q0, r0, P, Q, r and b.
    VECTOR xmami = xmax-xmin;
    VECTOR xmamieps = eeen*0.00001;
    xmami = max(xmami,xmamieps);
    VECTOR xmamiinv = eeen.pointdiv(xmami);
    VECTOR ux1 = upp-xval;
    VECTOR ux2 = ux1*ux1;
    VECTOR xl1 = xval-low;
    VECTOR xl2 = xl1*xl1;
    VECTOR uxinv = eeen.pointdiv(ux1);
    VECTOR xlinv = eeen.pointdiv(xl1);
    //
    
    VECTOR p0; p0 = zeron;
    VECTOR q0; q0 = zeron;
    p0 = max(df0,zeron);
    q0 = min(df0,zeron); q0*=-1;
    VECTOR pq0 = p0 + q0;
    p0 = p0 + pq0*0.001;
    q0 = q0 + pq0;
    p0 = p0 + xmamiinv*raa0;
    q0 = q0 + xmamiinv*raa0;
    p0 = p0*ux2;
    q0 = q0*xl2;
    prec r0 = f0 - p0.dot(uxinv) - q0.dot(xlinv);
    //
    int nCons = dg.nRow;
    std::vector<VECTOR> P(nCons);
    for (int i = 0; i < nCons; i++)
    {
        P[i] = max(dg[i], 0, n);
    }
    std::vector<VECTOR> Q(nCons);
    for (int i = 0; i < nCons; i++)
    {

        Q[i].initialize(n);
        Q[i] = min(dg[i], 0, n);
        Q[i] *= -1;
    }
    VECTOR b = VECTOR::zeros(nCons);
    for (int i = 0; i < nCons; i++)
    {
        VECTOR PQ = (P[i] + Q[i])*0.001 + xmamiinv*raa[i];
        
        P[i] += PQ; Q[i] += PQ;
        P[i] *= ux2;
        Q[i] *= xl2;

        b[i] += P[i].dot(uxinv) + Q[i].dot(xlinv);
    }
    b -= g;

    //
    // Solving the subproblem by a primal-dual Newton method
    subsolv(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,
            m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
    //
    // Calculations of f0app and fapp.
    ux1 = upp-xmma;
    xl1 = xmma-low;
    uxinv = eeen.pointdiv(ux1);
    xlinv = eeen.pointdiv(xl1);
    f0app = r0 + p0.dot(uxinv) + q0.dot(xlinv);
    for (int i = 0; i < nCons; i++)
    {

        fapp[i] = P[i].dot(uxinv) + Q[i].dot(xlinv) - b[i];
    }
    //
}

void OPTIMIZER::concheck(prec &conserv, int m,prec epsimin,prec &f0app,prec &f0valnew, VECTOR &fapp, VECTOR &fvalnew)
{
    conserv = 0;
    VECTOR eeem(m); eeem.reset(1.0);
    prec f0appe = f0app+epsimin;
    VECTOR fappe = fapp+eeem*epsimin;
    bool checked = true;
    if (f0appe >= f0valnew)
    {
        for (int i = 0; i < m; i++)
        {
            if (fappe[i] < fvalnew[i]) 
            {
                checked = false; 
                break;
            }
        }
        if (checked) conserv = 1;
    }
}

void OPTIMIZER::raaupdate(prec &raa0, VECTOR &raa,
                   VECTOR &xmma,VECTOR &xval,VECTOR &xmin,VECTOR &xmax,VECTOR &low,VECTOR &upp,prec f0valnew,VECTOR &fvalnew,prec f0app,VECTOR &fapp, prec raa0eps,
                   VECTOR &raaeps,prec epsimin)
{
    //
    int m = raa.length; int n = xmma.length;
    prec raacofmin = 1e-12;
    VECTOR eeem(m); eeem.reset(1.0);
    VECTOR eeen(n); eeen.reset(1.0);
    VECTOR xmami = xmax-xmin;
    VECTOR xmamieps = eeen* 0.00001;
    xmami = max(xmami,xmamieps);
    VECTOR xxux = (xmma-xval).pointdiv(upp-xmma);
    VECTOR xxxl = (xmma-xval).pointdiv(xmma-low);
    VECTOR xxul = xxux*xxxl;
    VECTOR ulxx = (upp-low).pointdiv(xmami);
    prec raacof = xxul.dot(ulxx);
    raacof = max(raacof,raacofmin);
    //
    prec f0appe = f0app + 0.5*epsimin;
    if (f0valnew > f0appe)
    {
        prec deltaraa0 = (1/raacof)*(f0valnew-f0app);
        prec zz0 = 1.1*(raa0 + deltaraa0);
        zz0 = min(zz0,10*raa0);
    //  zz0 = min(zz0,1000*raa0);
        raa0 = zz0;
    }
    //
    VECTOR fappe = fapp + eeem*(0.5*epsimin);
    VECTOR fdelta = fvalnew-fappe;
    VECTOR deltaraa = (fvalnew-fapp)*(1/raacof);
    VECTOR zzz = (raa + deltaraa)*1.1;
    zzz = min(zzz,raa*10);
    //zzz = min(zzz,1000*raa);
    for (int i = 0; i < m; i++)
    {
        if (fdelta[i] > 0) raa[i] = zzz[i];
    }
    //---------------------------------------------------------------------
}

//-----------------------------------------
// DEFINE INLET ELEMENTS AND AREA
//-----------------------------------------
void OPTIMIZER::define_inlet_pressure_elem()
{
    VECTOR_INT inlets;
    inlets = (*physics).inlet_bounds;

    int dim = (*physics).dim;
    int nBounds = (*physics).nBounds;
    
    MATRIX_INT inlet_elems((*physics).nElem, dim);
    int count_elem = 0;

    // BEGIN STREAMING
    std::ifstream bound_info_file_stream;
    bound_info_file_stream.open((*physics).bound_info_file_path);
    if (!bound_info_file_stream.is_open()) throw_line("ERROR, can't open input data file");
    std::string line;
    std::istringstream iss;

    getline(bound_info_file_stream, line); // skip general info

    for (int ibound = 0; ibound < nBounds; ibound++)
    {
        VECTOR_INT temp_info(2);
        STREAM::getRowVector(bound_info_file_stream, line, iss, temp_info);
        if (inlets.hasIn(temp_info[0]))
        {
            int temp_n_elem = temp_info[1];
            for (int iel = 0; iel < temp_n_elem; iel++)
            {
                VECTOR_INT tempElems(dim);
                STREAM::getRowVector(bound_info_file_stream, line, iss, tempElems);
                for (int inode = 0; inode < dim; inode++)
                {
                    inlet_elems[count_elem][inode] = tempElems[inode];
                    
                }
                count_elem += 1;
            }
        }
        else
        {
            STREAM::getLines(bound_info_file_stream, line, temp_info[1]);
        }
        getline(bound_info_file_stream, line);
    }  

    // CLOSE STREAMING
    bound_info_file_stream.close();

    inlet_elems.shrinkRows(count_elem); 
    (*physics).n_inlet_bounds_elems = count_elem;
    (*physics).inlet_bounds_elems = inlet_elems;
    //pause();

    // get Area
    VECTOR inlet_elems_surface(count_elem);

    MATRIX matCoord(dim,dim); 
    std::shared_ptr<prec*[]> coord = (*physics).coord.PP;
    for (int iel = 0; iel < count_elem; iel++)
    {
        for (int inode = 0; inode < dim; inode++)
        {
            int iglob = inlet_elems[iel][inode];
            for (int icomp = 0; icomp < dim; icomp++) matCoord[inode][icomp] = coord[iglob][icomp]; 
        }
        inlet_elems_surface[iel] = PHYSICS::get_surface(matCoord, dim);
    }
    (*physics).area_inlet_bounds_elems = inlet_elems_surface;


    // Get INLET NODES_V
    VECTOR_INT inlet_nodes_v;

    // BEGIN STREAMING
    std::ifstream bound_nodes_v_file_stream;
    bound_nodes_v_file_stream.open((*physics).bound_nodes_v_file_path);
    if (!bound_nodes_v_file_stream.is_open()) throw_line("ERROR, can't open input data file");
    // std::string line;
    // std::istringstream iss;

    getline(bound_nodes_v_file_stream, line); // skip general info

    for (int ibound = 0; ibound < nBounds; ibound++)
    {
        VECTOR_INT temp_info(2);
        STREAM::getRowVector(bound_nodes_v_file_stream, line, iss, temp_info);
        if (inlets.hasIn(temp_info[0]))
        {
            int temp_n_nodes_v = temp_info[1];
            VECTOR_INT tempNodes_v(temp_n_nodes_v);
            tempNodes_v.reset(0);
            STREAM::getColVector(bound_nodes_v_file_stream, line, iss, tempNodes_v, temp_n_nodes_v);
            inlet_nodes_v.append(tempNodes_v);
        }
        else
        {
            STREAM::getLines(bound_nodes_v_file_stream, line, temp_info[1]);
        }
        getline(bound_nodes_v_file_stream, line);
    }  
    
    // CLOSE STREAMING
    bound_nodes_v_file_stream.close();

    (*physics).inlet_nodes_v = inlet_nodes_v;
}

void OPTIMIZER::get_time_functional(int solution_time, MATRIX_INT &elem_v, VECTOR &Volume_v, MATRIX_INT &inlet_bound_elem, VECTOR &area_inlet_bound_elem)
{
    VECTOR temp_NS_sol = (*physics).NS_solution.get_row(solution_time);
    VECTOR temp_ADJ_sol = (*physics).NS_solution.get_row(solution_time);

    MATRIX U(dim, (*physics).nNodes_v);
    VECTOR P((*physics).nNodes);
    MATRIX Ua(dim, (*physics).nNodes_v);
    VECTOR Pa((*physics).nNodes);

    decompose_solution(temp_NS_sol, U, P);
    decompose_solution(temp_ADJ_sol, Ua, Pa);

    temp_obj_functional = 0.0;
    temp_func_val.resetZeros();
    temp_d_obj_functional.resetZeros();

    get_functional_and_opt_derivative(elem_v, Volume_v, U, Ua, inlet_bound_elem, area_inlet_bound_elem, P);       
}

void OPTIMIZER::get_functional_and_opt_derivative(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U, MATRIX &Ua,MATRIX_INT &elem, VECTOR &volume, VECTOR &P)
{
    eval_alpha_and_dAlpha();
    eval_functional(elem_v, volume_v, U, elem, volume, P);
    eval_opt_derivative(elem_v, volume_v, U, Ua);
}

void OPTIMIZER::eval_alpha_and_dAlpha()
{
    alpha.reset(alpha_min);
    prec alpha_factor = q*(alpha_max - alpha_min);
    prec dAlpha_factor = -q*(q+1)*(alpha_max-alpha_min);
    dAlpha.reset(dAlpha_factor / ((q+1)*(q+1)));
    for (int ioptnode = 0; ioptnode < n_nodes_in_dom; ioptnode++)
    {
        int iglob = nodeInDom[ioptnode];
        prec temp_gamma_acc = gamma_acc[ioptnode];
        prec temp_dgamma_acc = dgamma_acc[ioptnode];
        prec temp_gamma_plus_q = q + temp_gamma_acc;
        alpha[iglob] = alpha_min + alpha_factor * (1 - temp_gamma_acc) / temp_gamma_plus_q;
        dAlpha[iglob] = (dAlpha_factor / (temp_gamma_plus_q * temp_gamma_plus_q)) * temp_dgamma_acc;
    }
}

void OPTIMIZER::eval_functional(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U, MATRIX_INT &elem, VECTOR &area, VECTOR &P)
{
    temp_obj_functional = 0.0;
    temp_func_in_box = 0.0;
    eval_J_alpha(elem_v, volume_v, U);
    eval_J_gradU(elem_v, volume_v, U);
    eval_J_omega(elem_v, volume_v, U);
    eval_J_p_inlet(elem, area, P);
    temp_obj_functional = temp_func_val.sum();
}

void OPTIMIZER::eval_J_alpha(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U)
{
    temp_func_val[0] = 0.0;
    temp_no_weighted_func_val[0] = 0.0;

    if (abs(fWeights[0]) > 1e-12)
    {
        prec temp_alpha = 0.0;
        for (int iel = 0; iel < nElem_v; iel++)
        {
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[iel][iloc];
                prec tempFactor = 0.0;
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    prec Uval = U[icomp][iglob];
                    tempFactor += Uval*Uval;
                }
                prec temp_integral = volume_v[iel] / ((dim+1)*(dim+2)/2); 
                temp_alpha += alpha[iglob] * tempFactor * temp_integral; 
            }
        }
        temp_no_weighted_func_val[0] = temp_alpha;
        temp_func_val[0] = fWeights[0]*temp_alpha;

        // get func in top opt box
        prec temp_f_alpha = 0.0;
        for (int iel = 0; iel < n_elems_in_dom; iel++)
        {
            int globEl = elemInDom[iel];
            for (int iloc = 0; iloc < dim+1; iloc++)
            {
                int iglob = elem_v[globEl][iloc];
                // int optNode = optNodeFromGlobNode[iglob];
                prec tempFactor = 0;
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    prec Uval = U[icomp][iglob];
                    tempFactor += Uval*Uval;
                }
                prec temp_integral = volume_v[iel] / ((dim+1)*(dim+2)/2); 
                temp_f_alpha += alpha[iglob] * tempFactor * temp_integral;
            }
        }
        temp_func_in_box += temp_f_alpha;
    }
}

void OPTIMIZER::eval_J_gradU(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U)
{
    // grad
    temp_func_val[1] = 0.0;
    temp_no_weighted_func_val[1] = 0.0;

    prec temp_grad = 0.0;
    for (int iel = 0; iel < nElem_v; iel++)
    {
        MATRIX tempGrad;
        tempGrad.zeros(dim,dim);
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[jcomp][iel][iloc] + U[jcomp][iglob]*COEF[icomp][iel][iloc];
                }
            }
        }
        prec tempNorm = tempGrad.normFro();
        temp_grad += tempNorm*tempNorm * volume_v[iel];
    }
    // if (isnan(temp_grad)) 
    // {
    //     std::cout << "ERROR: temp_grad NAN\n";
    //     pause();
    // }
    temp_grad *= 0.5*mu;
    temp_no_weighted_func_val[1] = temp_grad;
    if (abs(fWeights[1]) > 1e-12) temp_func_val[1] = fWeights[1]*temp_grad;

    // get func in top opt box
    prec temp_f_grad = 0.0;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        MATRIX tempGrad;
        tempGrad.zeros(dim,dim);
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[jcomp][globEl][iloc] + U[jcomp][iglob]*COEF[icomp][globEl][iloc];
                }
            }
        }
        prec tempNorm = tempGrad.normFro();
        temp_f_grad += tempNorm*tempNorm * volume_v[globEl];
    }
    temp_f_grad *= 0.5*mu;
    if (abs(fWeights[1]) > 1e-12)
    {
        temp_f_grad *= fWeights[1];
    }
    else
    {
        temp_f_grad = 0.0;
    }
    temp_func_in_box += temp_f_grad;
}

void OPTIMIZER::eval_J_omega(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U) // eval the vorticity functional by calculating the vorticity
{
    //curl
    temp_func_val[2] = 0.0;
    temp_no_weighted_func_val[2] = 0.0;

    prec temp_omega = 0.0;
    VECTOR_INT ezComps;
    if (dim == 3)
    {
        ezComps.initialize(5);
        ezComps[0] = 0; ezComps[1] = 1; ezComps[2] = 2; ezComps[3] = 0; ezComps[4] = 1;
    }
    for (int iel = 0; iel < nElem_v; iel++)
    {
        VECTOR tempCurl(3);
        tempCurl.resetZeros();
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            if ( dim == 2)
            {
                tempCurl[0] = 0; tempCurl[1] = 0;
                tempCurl[2] += U[1][iglob] * COEF[0][iel][iloc] - U[0][iglob]*COEF[1][iel][iloc];
            }
            else if (dim == 3)
            {
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    int comp2 = ezComps[icomp+2];
                    int comp1 = ezComps[icomp+1];
                    tempCurl[icomp] += U[comp2][iglob] * COEF[comp1][iel][iloc] - U[comp1][iglob]*COEF[comp2][iel][iloc];
                }
            }
        }
        prec tempNorm = tempCurl.norm();
        temp_omega += tempNorm*tempNorm * volume_v[iel]; 
    }
    // if (isnan(temp_omega)) 
    // {
    //     std::cout << "ERROR: temp_omega NAN\n";
    //     pause();
    // }
    temp_omega *= 0.5*mu;
    temp_no_weighted_func_val[2] = temp_omega;
    if (abs(fWeights[2]) > 1e-12) temp_func_val[2] = fWeights[2]*temp_omega;

    // get func in top opt box
    prec temp_f_vort= 0.0;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        VECTOR tempCurl(3);
        tempCurl.resetZeros();
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            if ( dim == 2)
            {
                tempCurl[0] = 0; tempCurl[1] = 0;
                tempCurl[2] += U[1][iglob] * COEF[0][globEl][iloc] - U[0][iglob]*COEF[1][globEl][iloc];
            }
            else if (dim == 3)
            {
                for (int icomp = 0; icomp < dim; icomp++)
                {
                    int comp2 = ezComps[icomp+2];
                    int comp1 = ezComps[icomp+1];
                    tempCurl[icomp] += U[comp2][iglob] * COEF[comp1][globEl][iloc] - U[comp1][iglob]*COEF[comp2][globEl][iloc];
                }
            }
        }
        prec tempNorm = tempCurl.norm();
        temp_f_vort += tempNorm*tempNorm * volume_v[globEl]; 
    }
    temp_f_vort *= 0.5*mu;
    if (abs(fWeights[2]) > 1e-12)
    {
        temp_f_vort *= fWeights[2];
    }
    else
    {
        temp_f_vort = 0.0;
    }
    temp_func_in_box += temp_f_vort;
}

void OPTIMIZER::eval_J_omega_by_grad(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U) // eval the vorticity functional by calculating the anti-sym part of the gradient
{
    // curl
    temp_func_val[2] = 0.0;
    temp_no_weighted_func_val[2] = 0.0;

    prec temp_omega = 0.0;
    for (int iel = 0; iel < nElem_v; iel++)
    {
        MATRIX tempGrad;
        tempGrad.zeros(dim,dim);
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    tempGrad[icomp][jcomp] += U[icomp][iglob]*COEF[jcomp][iel][iloc] - U[jcomp][iglob]*COEF[icomp][iel][iloc];
                }
            }
            
        }
        prec tempNorm = tempGrad.normFro();
        temp_omega += tempNorm*tempNorm * volume_v[iel];
    }
    // if (isnan(temp_grad)) 
    // {
    //     std::cout << "ERROR: temp_grad NAN\n";
    //     pause();
    // }
    temp_omega *= 0.5 * 0.5*mu; //the first 0.5 is imposed because the norm of the vorticity is one half the norm of the anti-sym part of the gradient
    temp_no_weighted_func_val[2] = temp_omega;
    if (abs(fWeights[2]) > 1e-12) temp_func_val[2] = fWeights[2]*temp_omega;

    // get func in top opt box
    prec temp_f_omega = 0.0;
    for (int iel = 0; iel < n_elems_in_dom; iel++)
    {
        MATRIX tempGrad;
        tempGrad.zeros(dim,dim);
        int globEl = elemInDom[iel];
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[globEl][iloc];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                for (int jcomp = 0; jcomp < dim; jcomp++)
                {
                    tempGrad[icomp][jcomp] += U[icomp][iglob] * COEF[jcomp][globEl][iloc] - U[jcomp][iglob]*COEF[icomp][globEl][iloc];
                }
            }
        }
        prec tempNorm = tempGrad.normFro();
        temp_f_omega += tempNorm*tempNorm * volume_v[globEl];
    }
    temp_f_omega *= 0.5 * 0.5*mu;
    if (abs(fWeights[2]) > 1e-12)
    {
        temp_f_omega *= fWeights[2];
    }
    else
    {
        temp_f_omega = 0.0;
    }
    temp_func_in_box += temp_f_omega;
}

void OPTIMIZER::eval_J_p_inlet(MATRIX_INT &elem, VECTOR &area, VECTOR &P)
{
    // inlet pressure
    temp_func_val[3] = 0.0;
    temp_no_weighted_func_val[3] = 0.0;

    int n_nodes_x_bound_el = (*physics).dim;

    int n_elems = (*physics).n_inlet_bounds_elems;

    prec temp_pressure = 0.0;
    for (int iel = 0; iel < n_elems; iel++)
    {
        prec temp_el_pressure = 0.0;
        for (int iloc = 0; iloc < dim; iloc++)
        {
            int iglob = elem[iel][iloc];
            temp_el_pressure += P[iglob];
        }
        temp_el_pressure *= area[iel] / n_nodes_x_bound_el;
        temp_pressure += temp_el_pressure;
    }
    // if (isnan(temp_pressure)) 
    // {
    //     std::cout << "ERROR: temp_pressure NAN\n";
    //     pause();
    // }
    temp_no_weighted_func_val[3] = temp_pressure;
    if (abs(fWeights[3]) > 1e-12) temp_func_val[3] = fWeights[3]*temp_pressure;
    
    // get func in top opt box
    prec temp_f_p = 0.0;
    for (int iel = 0; iel < n_elems; iel++)
    {
        prec temp_el_pressure = 0.0;
        for (int iloc = 0; iloc < dim; iloc++)
        {
            int iglob = elem[iel][iloc];
            if (nodeInDom.hasIn(iglob))
            {
                temp_el_pressure += P[iglob];
            }                
        }
        temp_el_pressure *= area[iel] / n_nodes_x_bound_el;
        temp_f_p += temp_el_pressure;
    }
    if (abs(fWeights[3]) > 1e-12)
    {
        temp_f_p *= fWeights[3];
    }
    else
    {
        temp_f_p = 0.0;
    }
    temp_func_in_box += temp_f_p;
}

void OPTIMIZER::eval_opt_derivative(MATRIX_INT &elem_v, VECTOR &volume_v, MATRIX &U, MATRIX &Ua)
{
    temp_d_obj_functional.resetZeros();

    VECTOR tempU(dim);
    VECTOR tempUa(dim);

    for (int iel = 0; iel < nElem_v; iel++)
    {
        for (int iloc = 0; iloc < dim+1; iloc++)
        {
            int iglob = elem_v[iel][iloc];
            for (int icomp = 0; icomp < dim; icomp++)
            {
                tempU[icomp] = U[icomp][iglob];
                tempUa[icomp] = Ua[icomp][iglob];
            }
            temp_d_obj_functional[iglob] += (dAlpha[iglob] * tempU.dot(tempUa)) / (dim+1) * volume_v[iel]; // * alpha?
            temp_d_obj_functional[iglob] += fWeights[0]*(dAlpha[iglob] * tempU.dot(tempU)) / (dim+1) * volume_v[iel];
        }
    }
}

void OPTIMIZER::eval_gamma_acc_and_derivative()
{
    gamma_acc.resetZeros();
    dgamma_acc.resetZeros();
    switch (gamma_acc_case)
    {
        case 0: // no projection
        {   
            for (int inode = 0; inode < n_nodes_in_dom; inode++)
            {
                gamma_acc[inode] = gamma_filter[inode];
                dgamma_acc[inode] = dgamma_filter[inode];
            }
            break;
        }
        case 1: //projection function: -2x^3 + 3x^2
        {
            for (int inode = 0; inode < n_nodes_in_dom; inode++)
            {
                prec gamma_filter_node = gamma_filter[inode];
                prec dgamma_filter_node = dgamma_filter[inode];
                gamma_acc[inode] = -2*(gamma_filter_node*gamma_filter_node*gamma_filter_node) + 3*gamma_filter_node*gamma_filter_node;
                dgamma_acc[inode] = 6*gamma_filter_node*(1 - gamma_filter_node) * dgamma_filter_node;
            }
            break;
        }
        case 2: // projection function: 0.5 + tanh(b*(x-0.5))/tanh(b*0.5)
        {
            prec temp_change = (*physics).gamma_change;
            eval_beta_for_projection_filter(temp_change);
            gamma_acc_mean_value = (*physics).gamma_max / 2;
            prec gamma_acc_mean_value_min = 0.2;
            if ((gamma_acc_mean_value < 0.0) || (gamma_acc_mean_value > 0.5))
            {
                throw_line("ERROR: invalid mean value of the gamma projector.\n");
            }
            else if (gamma_acc_mean_value < gamma_acc_mean_value_min)
            {
                gamma_acc_mean_value = gamma_acc_mean_value_min;
            }
            prec factor = tanh(beta_proj*gamma_acc_mean_value) + tanh(beta_proj*(1-gamma_acc_mean_value));
            prec constant = tanh(beta_proj*gamma_acc_mean_value);
            for (int inode = 0; inode < n_nodes_in_dom; inode++)
            {
                prec gamma_filter_node = gamma_filter[inode];
                prec dgamma_filter_node = dgamma_filter[inode];
                prec temp_arg = beta_proj*(gamma_filter_node-gamma_acc_mean_value);
                gamma_acc[inode] = (constant + tanh(temp_arg)) / factor; // [tanh(beta*alpha) + tanh*beta*(gamma_node-alpha)] / [tanh(beta*alpha) + tanh(beta*(1-alpha))]
                dgamma_acc[inode] = (beta_proj / factor) * (1.0 / (cosh(temp_arg)*cosh(temp_arg))) * dgamma_filter_node; // (beta / (2*tanh(0.5*beta))) * (1 / (cosh(beta*(gamma_node-0.5))^2)) * dgamma_filter/dgamma
            }
            break;
        }
        default:
            {
                throw_line("ERROR: not handled gamma acceleration case\n");
                break;
            }
    }
}

void OPTIMIZER::eval_beta_for_projection_filter(prec temp_change)
{
    // EVAL q_beta
    if (beta_interpolation == 0)
    {
        q_beta = 0;
    }
    else
    {
        if ((*physics).turn_on_gamma_acc == 0)
        {
            q_beta = (crit_change - change_min) * (crit_beta - beta_min_init) / ((beta_max-beta_min_init)*(change_max-crit_change) - (crit_beta-beta_min_init)*(change_max-change_min)); // defined to have beta(cirt_change) = crit_beta;
            if (q_beta <= 0)
            {
                throw_line("ERROR: beta slope coeffcient in gamma acceleration has negative value.\n");
            }
        }
    }

    // Gamma accelerator activation
    if (((*physics).turn_on_gamma_acc == 0) && ((temp_change > change_max) || (*physics).vol_fract <= 0.8))
    {
        std::cout << "\n\n## Optimization Log: | Gamma projection turned ON |\n\n";
        (*physics).turn_on_gamma_acc = 1;
    }

    // Eval beta for gamma projection filter
    if ((*physics).turn_on_gamma_acc == 0)
    {
        beta_proj = beta_min; // almost linear projection <=> no projection <=> no gamma acceleration
    }
    else if ((*physics).turn_on_gamma_acc >= 1)
    {
        if (temp_change > change_max)
        {
            beta_proj = beta_min;
            (*physics).turn_on_gamma_acc = 1;
            //std::cout << "\n\nq: "<< q_beta << "\tBETA: " << beta << "\tCASE: " << (*physics).turn_on_gamma_acc << "\n\n";
        }
        else if (temp_change < change_min)  
        {
            beta_proj = beta_max;
            (*physics).turn_on_gamma_acc = 3;
            //std::cout << "\n\nq: "<< q_beta << "\tBETA: " << beta << "\tCASE: " << (*physics).turn_on_gamma_acc << "\n\n";
        }
        else
        {
            beta_proj = beta_min + (beta_max - beta_min) * q_beta * (change_max - temp_change) / ((change_max-change_min)*q_beta + temp_change - change_min);
            (*physics).turn_on_gamma_acc = 2;
            //std::cout << "\n\nq: "<< q_beta << "\tBETA: " << beta << "\tCASE: " << (*physics).turn_on_gamma_acc << "\n\n";
        }
    }
    else
    {
        throw_line("ERROR: non handled turn_on_gamma_acc value.\n");
    }

    // Update beta_min_threshold or beta_proj to guarentee smooth convergence
    if (beta_proj > beta_min_threshold)
    {
        beta_min_threshold = beta_proj;
    }
    else
    {
        beta_proj = beta_min_threshold;
    }
    // Update beta_min_threshold or beta_proj to guarentee smooth convergence
    // if (beta_proj > beta_min)
    // {
    //     beta_min = beta_proj;
    // }
    
}

void OPTIMIZER::topology_optimization_diffusive_filter()
{
    gamma_filter.resetZeros();
    dgamma_filter.resetZeros();
    switch (diffusion_filter_case)
    {
        case 0:
        {
            for (int inode = 0; inode < n_nodes_in_dom; inode++)
            {
                gamma_filter[inode] = gamma[inode];
                dgamma_filter[inode] = 1.0;
            }
            //std::cout << "filter case: 0\n";
            break;
        }
        case 1: case 2: // filters based on pythagorean means (1: artihmetic, 2:geometric)
        {
            diffusionFilter.filter_gamma_and_derivative(gamma, gamma_filter, dgamma_filter);
            break;
        }
        default:
        {
            throw_line("ERROR: non valid diffusive filter for the top. opt. parameter case\n");
            break;
        }
    }
}



