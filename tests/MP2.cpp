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


MP2::MP2() : HartreeFock(0.01, 0.01){

    READIN::TEI("data/eri.dat", &TEI_MO);
    READIN::TEI("data/eri.dat", &TEI_AO);

    HF->Iterate();

}

void MP2::NoddyAlgorithm () {
    /**
     AO to MO integral transformation using a single N^8 step
     Both Two Electron Integrals are stored in arrays, taking advantage of
     Permuational Symmetry
     **/

    int i, j, k, l, ijkl;
    int p, q, r, s, pq, rs, pqrs;

    for(i=0, ijkl=0; i < NUM_ORB; i++) { // instantiate compound and i
    for(j=0; j <= i; j++) { // iterate over unique sets of j (less than i)
    for(k=0; k <= i; k++) { // iterate over unique sets of k (less than i)
    for(l=0; l <= (i==k ? j : k); l++, ijkl++) { // iterate over unique l

    for(p=0; p < NUM_ORB; p++) { // Over all orbitals
    for(q=0; q < NUM_ORB; q++) { // Over all orbitals

    pq = INDEX(p,q);

        for(r=0; r < NUM_ORB; r++) { // Over all orbitals
        for(s=0; s < NUM_ORB; s++) { // Over all orbitals

        rs = INDEX(r,s);
        pqrs = INDEX(pq,rs);
        TEI_MO(ijkl) += HF->C0(p,i) * HF->C0(q,j) * HF->C0(r,k) * HF->C0(s,l) * TEI_AO(pqrs);

        }}
    }}}}}}
  }
