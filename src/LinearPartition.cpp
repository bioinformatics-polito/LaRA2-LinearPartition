/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h> 

#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "bpp.cpp"

#define SPECIAL_HP

using namespace std;

namespace linearp
{
template <typename pf_type, typename value_type>
unsigned long BeamCKYParser<pf_type, value_type>::quickselect_partition(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper) {
    pf_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
template <typename pf_type, typename value_type>
pf_type BeamCKYParser<pf_type, value_type>::quickselect(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

template <typename pf_type, typename value_type>
pf_type BeamCKYParser<pf_type, value_type>::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        pf_type newalpha = (k >= 0 ? bestC[k].alpha : pf_type(0.0)) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if ((int)scores.size() <= beam) return VALUE_MIN;
    pf_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

template <typename pf_type, typename value_type>
void BeamCKYParser<pf_type, value_type>::prepare(unsigned len) {
    seq_length = len;

    nucs = new int[seq_length];
    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<int, State>[seq_length];
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
    
    scores.reserve(seq_length);
}

template <typename pf_type, typename value_type>
void BeamCKYParser<pf_type, value_type>::postprocess() {

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
}

template <typename pf_type, typename value_type>
void BeamCKYParser<pf_type, value_type>::parse(string& seq) {
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (unsigned i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    vector<int> next_pair[NOTON];
    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            // next_pair
            next_pair[nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = seq_length-1; j >=0; --j) {
                next_pair[nuci][j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }
        }
    }

#ifdef SPECIAL_HP
    if (useVienna)
        v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif

    if (useVienna) {
        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;
    } else {
        if(seq_length > 0) Fast_LogPlusEquals(bestC[0].alpha, score_external_unpaired(0, 0));
        if(seq_length > 1) Fast_LogPlusEquals(bestC[1].alpha, score_external_unpaired(0, 1));
    }

    value_type newscore;
    for(unsigned j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && (int)beamstepH.size() > beam) beam_prune(beamstepH);

            {
                // for nucj put H(j, j_next) into H[j_next]
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;
                    if (useVienna) {
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
                    } else {
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore);
                    }
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
                for (auto &item : beamstepH) {
                    unsigned i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)=
                        if (useVienna) {
                            int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                            if (jnext-i-1 == 4) // 6:tetra
                                tetra_hex_tri = if_tetraloops[i];
                            else if (jnext-i-1 == 6) // 8:hexa
                                tetra_hex_tri = if_hexaloops[i];
                            else if (jnext-i-1 == 3) // 5:tri
                                tetra_hex_tri = if_triloops[i];
#endif
                            newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                            Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
                        } else {
                            newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
                            Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore);
                        }
                    }

                    // 2. generate p(i, j)
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && (int)beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
                        if (useVienna)
                            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
                        else {
                            newscore = score_multi_unpaired(j, jnext - 1);
                            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha + newscore);
                        }
                    }
                }

                // 2. generate P (i, j)
                {
                    if (useVienna) {
                        newscore = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                        Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore/kT);
                    } else {
                        newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                        Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore);
                    }
                }
            }
        }

        // beam of P
        {   
            if (beam > 0 && (int)beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    value_type precomputed = 0;
                    if (!useVienna)
                        precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
                        
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && (unsigned)q == j + 1) {
                                // helix
                                if (useVienna) {
                                    newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                                 nuci_1, nuci, nucj, nucj1);
    
                                    // SHAPE for Vienna only
                                    if (use_shape)
                                    {
                                        newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                    }
    
    
                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore/kT);
                                } else {
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);
                                }
                            } else {
                                // single branch
                                if (useVienna) {
                                    newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                                nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore/kT);
                                } else {
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                               precomputed +
                                               score_single_without_junctionB(p, q, i, j,
                                                                              nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);
                                }
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                    if (useVienna) {
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore/kT);
                    } else {
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore);
                    }
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
                    pf_type m1_alpha;
                    if (useVienna) {
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        m1_alpha = state.alpha + newscore/kT;
                    } else {
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        m1_alpha = state.alpha + newscore;
                    }
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        State& prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        if (useVienna) {
                            newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length);
                            Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore/kT);      
                        } else {
                            newscore = score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                            Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore);
                        }
                    } else {
                        if (useVienna) {
                            newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                            Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore/kT);       
                        } else {
                            newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                            Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore);
                        }
                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && (int)beamstepM2.size() > beam) beam_prune(beamstepM2);

            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                    int nucp = nucs[p];
                    int q = next_pair[nucp][j];
                    if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
                        if (useVienna)
                            Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha);      

                        else {
                            newscore = score_multi_unpaired(p+1, i-1) +
                                            score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha + newscore);      
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
            }
        }
        
        // beam of M
        {
            if (beam > 0 && (int)beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    if (useVienna)
                        Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
                        
                    else {
                        newscore = score_multi_unpaired(j + 1, j + 1);
                        Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha + newscore); 
                    }
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                if (useVienna)
                    Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
                    
                else {
                    newscore = score_external_unpaired(j+1, j+1);
                    Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha + newscore); 
                }
            }
        }
    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    // unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;

    if (useVienna)
        m_energy = -kT * viterbi.alpha / 100.0; // Free Energy of Ensemble kcal/mol
    else
        m_logPartitionCoef = viterbi.alpha;

    if(is_verbose) fprintf(stderr,"Partition Function Calculation Time: %.2f seconds.\n", parse_elapsed_time);

    fflush(stdout);

    // lhuang
    if(pf_only && !forest_file.empty()) dump_forest(seq, true); // inside-only forest

    if(!pf_only){
        outside(next_pair);
    	if (!forest_file.empty())
    	  dump_forest(seq, false); // inside-outside forest
            cal_PairProb(viterbi);

        if (mea_) PairProb_MEA(seq);

        if (threshknot_) ThreshKnot(seq);
    }
    postprocess();
    return;
}

template <typename pf_type, typename value_type>
void BeamCKYParser<pf_type, value_type>::print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        State & state = item.second;
        if (inside_only) fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
            fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

template <typename pf_type, typename value_type>
void BeamCKYParser<pf_type, value_type>::dump_forest(string seq, bool inside_only) {  
    printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
    FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
    fprintf(fptr, "%s\n", seq.c_str());
    int n = seq.length(), j;
    for (j = 0; j < n; j++) {
        if (inside_only) fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
        else fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
    }
    double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
    for (j = 0; j < n; j++) 
        print_states(fptr, bestP[j], j, "P", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM[j], j, "M", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

template <typename pf_type, typename value_type>
BeamCKYParser<pf_type, value_type>::BeamCKYParser(bool useVienna,
                                                  int beam_size,
                                                  bool nosharpturn,
                                                  float bppcutoff,
                                                  bool verbose,
                                                  string bppfile,
                                                  string bppfileindex,
                                                  bool pfonly,
                                                  string forestfile,
                                                  bool mea,
                                                  float MEA_gamma,
                                                  string MEA_file_index,
                                                  bool MEA_bpseq,
                                                  bool ThreshKnot,
                                                  float ThreshKnot_threshold,
                                                  string ThreshKnot_file_index,
                                                  string shape_file_path)
    : useVienna(useVienna),
      beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      bpp_cutoff(bppcutoff),
      is_verbose(verbose),
      bpp_file(bppfile),
      bpp_file_index(bppfileindex),
      pf_only(pfonly),
      forest_file(forestfile), 
      mea_(mea),
      gamma(MEA_gamma),
      mea_file_index(MEA_file_index),
      bpseq(MEA_bpseq),
      threshknot_(ThreshKnot),
      threshknot_threshold(ThreshKnot_threshold),
      threshknot_file_index(ThreshKnot_file_index){
    if(useVienna)
        initialize();
    else {
        initialize();
        initialize_cachesingle();
    }

    if (shape_file_path != "" ){
        use_shape = true;
        int position;
        string data;

        double temp_after_mb_shape;

        ifstream in(shape_file_path);

        if (!in.good()){
            cout<<"Reading SHAPE file error!"<<endl;
            assert(false);
        }

        // actually, we can combine the SHAPE_data and the energy_stack together
        while (!(in >> position >> data).fail()) {
            // cout<<"position data "<< int(position)<<endl<<data<<endl;
            // assert(int(position) == SHAPE_data.size() + 1);
            // cout<<"data "<<data<<endl;
            if (isdigit(int(data[0])) == 0){
                SHAPE_data.push_back(double((-1.000000)));
            }

            else {
                SHAPE_data.push_back(stod(data));
            }
            

        }

        for (size_t i = 0; i<SHAPE_data.size(); i++){
            temp_after_mb_shape = SHAPE_data[i] < 0 ? 0. : (m * log(SHAPE_data[i] + 1) + b);

            pseudo_energy_stack.push_back((int)roundf(temp_after_mb_shape * 100.));

            assert(pseudo_energy_stack.size() == i + 1 );

            // cout<<"pseudo energy "<<i<<' '<<SHAPE_data[i]<<' '<<temp_after_mb_shape<<' '<<pseudo_energy_stack[i]<<' '<<pseudo_energy_stack.size()<<endl;

        }
    }
}

template class BeamCKYParser<viennaModel::pf_type, viennaModel::value_type>;
template class BeamCKYParser<contraModel::pf_type, contraModel::value_type>;

} // end namespace linearp
