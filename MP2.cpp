//
//  MP2.cpp
//  Quantum
//
//  Created by Isabel Craig on 1/26/18.
//  Copyright © 2018 Isabel Craig. All rights reserved.
//

#include "MP2.hpp"
#include "HartreeFock.hpp"
#include "Read.hpp"


void NoddyAlgorithm (MulikenMatrix *TEI_MO, MulikenMatrix *TEI_AO) {
    /**
        AO to MO integral transformation using a single N^8 step
        Both Two Electron Integrals are stored in arrays, taking advantage of
        Permuational Symmetry
     **/
    
    int i, j, k, l, ijkl;
    int p, q, r, s, pq, rs, pqrs;
    
    for(i=0, ijkl=0; i < nao; i++) { // instantiate compound and i
    for(j=0; j <= i; j++) { // iterate over unique sets of j (less than i)
    for(k=0; k <= i; k++) { // iterate over unique sets of k (less than i)
    for(l=0; l <= (i==k ? j : k); l++, ijkl++) { // iterate over unique l
    
    for(p=0; p < nao; p++) {
    for(q=0; q < nao; q++) {
        
        pq = INDEX(p,q);
        
        for(r=0; r < NUM_ORB; r++) {
        for(s=0; s < NUM_ORB; s++) {
            
            rs = INDEX(r,s);
            pqrs = INDEX(pq,rs);
            (*TEI_MO)(ijkl) += C[p][i] * C[q][j] * C[r][k] * C[s][l] * (*TEI_AO)(pqrs);
            
        }}
    }}}}}}
}
