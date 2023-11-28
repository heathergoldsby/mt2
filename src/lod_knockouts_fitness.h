
#ifndef _EALIFE_LOD_KNOCKOUTS_FITNESS_H_
#define _EALIFE_LOD_KNOCKOUTS_FITNESS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
//#include <ea/analysis/tool.h>
#include <ea/analysis.h>
//#include <ea/functional.h>
#include <ea/digital_evolution/instruction_set.h>
//#include <ea/digital_evolution/discrete_spatial_environment.h>
#include <ea/digital_evolution/environment.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/math/information.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "gls.h"

LIBEA_MD_DECL(LOD_START_ANALYSIS, "ea.mt.lod_start_analysis", int);
LIBEA_MD_DECL(LOD_END_ANALYSIS, "ea.mt.lod_end_analysis", int);
LIBEA_MD_DECL(ANALYSIS_LOD_REPS, "ea.mt.lod_analysis_reps", int);
LIBEA_MD_DECL(ANALYSIS_LOD_START_COST, "ea.mt.lod_start_cost", int);
LIBEA_MD_DECL(ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE, "ea.mt.lod_timepoint_to_analyze", int);
LIBEA_MD_DECL(ANALYSIS_LOD_START_REP, "ea.mt.lod_analysis_start_rep", int);
LIBEA_MD_DECL(ANALYSIS_MUTATIONS_OFF, "ea.mt.lod_analysis_mutations_off", int);
LIBEA_MD_DECL(ONLY_MC, "ea.mt.only_mc", int);


namespace ealib {
    namespace analysis {
        
        
    
    LIBEA_ANALYSIS_TOOL(lod_fitness_combo) {
        
        datafile df("lod_fitness.dat");
        df.add_field("timepoint")
        .add_field("mc_or_uni")
        .add_field("count")
        .add_field("iteration")
        .add_field("time_to_fill")
        .add_field("workload")
        .add_field("num_org")
        .add_field("total_cells")
        ;
        
        
        datafile df2("lod_fit_summary.dat");
        df2.add_field("timepoint")
        .add_field("num_unicell_revertants")
        .add_field("num_viable_unicells")
        .add_field("num_inviable_unicells")
        .add_field("update")
        ;
        
        
        datafile df3("multicell_detail.dat");
        df3.add_field("timepoint")
        .add_field("count")
        .add_field("iteration")
        .add_field("multicell")
        .add_field("organism_id")
        .add_field("organism_parent_id")
        .add_field("cell")
        .add_field("workload")
        .add_field("propagule_eligible")
        .add_field("cell_resources")
        .add_field("not")
        .add_field("nand")
        .add_field("and")
        .add_field("ornot")
        .add_field("or")
        .add_field("andnot")
        .add_field("nor")
        .add_field("xor")
        .add_field("equals")
        .add_field("update")
        .add_field("multicell_resources")
        .add_field("total_workload")
        .add_field("genome")
        ;
        
        
                             
        datafile df4("unicell_detail.dat");
        df4.add_field("timepoint")
        .add_field("genome_location")
        .add_field("mutation")
        .add_field("iteration")
        .add_field("multicell")
        .add_field("organism_id")
        .add_field("organism_parent_id")
        .add_field("cell")
        .add_field("workload")
        .add_field("propagule_eligible")
        .add_field("cell_resources")
        .add_field("not")
        .add_field("nand")
        .add_field("and")
        .add_field("ornot")
        .add_field("or")
        .add_field("andnot")
        .add_field("nor")
        .add_field("xor")
        .add_field("equals")
        .add_field("update")
        .add_field("multicell_resources")
        .add_field("total_workload")
        .add_field("genome")
        ;
        
        datafile df5("inviable_unicell_detail.dat");
        df5.add_field("timepoint")
        .add_field("count")
        .add_field("iteration")
        .add_field("workload")
        .add_field("propagule_eligible")
        .add_field("cell_resources")
        .add_field("not")
        .add_field("nand")
        .add_field("and")
        .add_field("ornot")
        .add_field("or")
        .add_field("andnot")
        .add_field("nor")
        .add_field("xor")
        .add_field("equals")
        .add_field("update")
        .add_field("multicell_resources")
        .add_field("total_workload")
        ;

        
        int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);
        int track_details = get<TRACK_DETAILS>(ea,1);
        int num_rep = get<ANALYSIS_LOD_REPS>(ea,1);
        int mutations_off = get<ANALYSIS_MUTATIONS_OFF>(ea,0);
        int only_mc = get<ONLY_MC>(ea,0);
        line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
        //        typename line_of_descent<EA>::iterator i=lod.begin(); i++;
        
        
        
        
         
        typename line_of_descent<EA>::iterator i;
        int lod_step = 0;
        
        if (timepoint == 1) {
            i = lod.end(); --i;
        } else if (timepoint == 0) {
            i=lod.begin(); i++;
            lod_step++;
            // find the first to transition
            for( ; i!=lod.end(); i++) {
                lod_step++;
                if ((timepoint != 0) && (timepoint == lod_step)) {
                    break;
                }
                if ((timepoint == 0) && (i->size() > 2)) {
                    break;
                }
            }
        } else {
            i=lod.begin(); i++;
            lod_step++;
            for( ; i!=lod.end(); i++) {
                lod_step++;
                if (timepoint == lod_step) {
                    break;
                }
            }
        }
                
        
        
        // for 100 reps we run the multicell and gather info. 
        for (int nr = 0; nr < num_rep; nr++) {
            
            // should define checkpoint + analysis input
            //ea is the thing loaded from the checkpoint; EA is its type
            EA metapop; // a new EA
            //typedef metadata md_type;
            
            typename EA::md_type md(ea.md());
            // override md settings (pop size, geo, etc)
            
            
            metapop.initialize(md);
            put<METAPOPULATION_SIZE>(100, metapop);
            put<RUN_UPDATES>(10000, metapop);
            int new_seed = ea.rng().uniform_integer();
            put<RNG_SEED>(new_seed, metapop);
            metapop.reset_rng(new_seed);
            
            if (mutations_off) {
                put<TASK_MUTATION_PER_SITE_P>(0, metapop);
                //put<MUTATION_PER_SITE_P>(0, metapop);
                //put<GERM_MUTATION_PER_SITE_P>(0, metapop);
            }
        
            
            typename EA::individual_ptr_type control_mc = ea.make_individual(*i->traits().founder());
            control_mc->traits()._founder = i->traits().founder();
            
            put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *control_mc);
            
            if (mutations_off) {
                put<TASK_MUTATION_PER_SITE_P>(0, *control_mc);
                //put<MUTATION_PER_SITE_P>(0, *control_mc);
                //put<GERM_MUTATION_PER_SITE_P>(0, *control_mc);
            }
            
            new_seed = ea.rng().uniform_integer();
            put<RNG_SEED>(new_seed, *control_mc);
            control_mc->reset_rng(new_seed);
                        
            typename EA::population_type init_mc;
            init_mc.insert(init_mc.end(),control_mc);
            
            std::swap(metapop.population(), init_mc);
            
            add_event<mt_gls_propagule>(metapop);
            add_event<subpopulation_founder_event>(metapop);

            int max_size = 32;
            int max_update = 50000;
            int cur_update = 0;
            
            while ((metapop.size() < max_size) &&
                   (cur_update < max_update)){
                metapop.update();
                ++cur_update;
            }
            
            if (track_details){
            std::string new_filename = "mt_gls_detail_mc_" + std::to_string(nr) + ".dat";
            std::rename("mt_gls_detail.dat", new_filename.c_str());
            }
                    
            // get workload
            float total_workload = 0;
            float total_cells = 0;
            int multicell_count = 0;
            typedef typename EA::subpopulation_type::population_type subpop_type;
            
            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                total_cells += j->size();
                int cell_count = 0;
                int org_id = get<ORGANISM_ID>(*j,0);
                int org_parent_id = get<ORGANISM_PARENT_ID>(*j, 0);
                
            
                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                    
                    typename EA::subpopulation_type::individual_type& org=**m;
                    total_workload += get<WORKLOAD>(org, 0.0);
                    if (track_details) {
                        
                        /*         datafile df3("multicell_detail.dat");
                        df3.add_field("timepoint")
                        .add_field("count")
                        .add_field("iteration")
                        .add_field("multicell")
                        .add_field("organism_id")
                        .add_field("organism_parent_id")
                        .add_field("cell")
                        .add_field("workload")
                        .add_field("propagule_eligible")
                        .add_field("cell_resources")
                        .add_field("not")
                        .add_field("nand")
                        .add_field("and")
                        .add_field("ornot")
                        .add_field("or")
                        .add_field("andnot")
                        .add_field("nor")
                        .add_field("xor")
                        .add_field("equals")
                        .add_field("update")
                        .add_field("multicell_resources")
                        .add_field("total_workload")
                        .add_field("genome")
                        ;
                         */
                        df3.write(timepoint)
                        .write(0)
                        .write(nr)
                        .write(multicell_count)
                        .write(org_id)
                        .write(org_parent_id)
                        .write(cell_count)
                        .write(get<WORKLOAD>(org, 0.0))
                        .write(get<GERM_STATUS>(org,true))
                        .write(get<SAVED_RESOURCES>(org,0.0))
                        .write(get<TASK_NOT>(org,0))
                        .write(get<TASK_NAND>(org,0))
                        .write(get<TASK_AND>(org,0))
                        .write(get<TASK_ORNOT>(org,0))
                        .write(get<TASK_OR>(org,0))
                        .write(get<TASK_ANDNOT>(org,0))
                        .write(get<TASK_NOR>(org,0))
                        .write(get<TASK_XOR>(org,0))
                        .write(get<TASK_EQUALS>(org,0))
                        .write(cur_update)
                        .write(get<GROUP_RESOURCE_UNITS>(*j,0.0))
                        .write(total_workload);
                        
                        df3.write("\"");
                        for(typename EA::subpopulation_type::genome_type::iterator k2=org.genome().begin(); k2!=org.genome().end(); ++k2) {
                            df3.write(*k2)
                            .write(" ");
                        }
                        df3.write("\"");
                        df3.endl();
                        
                        
                    }
                    cell_count++;
                }
                multicell_count++;
            }
            
            
            df.write(timepoint)
            .write("mc")
            .write("0")
            .write(nr)
            .write(cur_update)
            .write(total_workload)
            .write(metapop.size())
            .write(total_cells);
            df.endl();
        }
        
        if (only_mc) {
            return;
        }
        
        
        // **i is the EA, AS OF THE TIME THAT IT DIED!
        
        // To replay, need to create new eas for each knockout exper.
        // setup the population (really, an ea):
        
        typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
        int birth_up = get<IND_BIRTH_UPDATE>(*i->traits().founder(),0);
        
        // ok we need to iterate through size...
        // fixed size 100 genome...
        int uni_count = 0;
        int num_uni_viable = 0;
        int num_uni_inviable = 0;
        int num_uni = 0;
        

        // for all places in the genome, we ...
        for (int z =0; z < knockout_loc->population()[0]->genome().size(); z++) {
            // iterate through each instruction in the isa
            for (int q = 0; q < control_ea->isa().size(); q++) {
                typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                
                put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                
                
                knockout_loc->population()[0]->genome()[z] = q;
                
                // we set our mutation rates to 0 here to see what happens in the absence of mutations.
                put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc);
                put<MUTATION_PER_SITE_P>(0, *knockout_loc);
                put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc);
                int new_seed = ea.rng().uniform_integer();
                put<RNG_SEED>(new_seed, *knockout_loc);
                knockout_loc->reset_rng(new_seed);
                
                int cur_update = 0;
                int update_max = 10000;
                // and run till the group amasses the right amount of resources
                while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                       (cur_update < update_max)){
                    knockout_loc->update();
                    ++cur_update;
                }
                
                
                // if it is a unicell in the absence of mutations, we...
                if (knockout_loc->population().size() < 2) {
                    num_uni++;
                    float total_workload = 0;
                    if (cur_update == update_max) {
                        num_uni_inviable++;
                        continue;
                    }
                    
                    
                    num_uni_viable++;
                    
                    // for each unicell, we recreate it. Then we knock it out
                    for (int nr = 0; nr < num_rep; nr++) {
                        typename EA::individual_ptr_type knockout_loc2 = ea.make_individual(*i->traits().founder());

                        put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc2);
                        put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc2);
               
                        
                        knockout_loc2->population()[0]->genome()[z] = q;
                        
                        if (mutations_off) {
                            put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc2);
                            //put<MUTATION_PER_SITE_P>(0, *knockout_loc2);
                            //put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc2);
                        }
                        
                        
                        // should define checkpoint + analysis input
                        //ea is the thing loaded from the checkpoint; EA is its type
                        EA metapop; // a new EA
                        //typedef metadata md_type;
                        typename EA::md_type md(ea.md());
                        // override md settings (pop size, geo, etc)
                        
                        
                        metapop.initialize(md);
                        put<METAPOPULATION_SIZE>(32, metapop);
                        put<RUN_UPDATES>(10000, metapop);
                        new_seed = ea.rng().uniform_integer();
                        put<RNG_SEED>(new_seed, metapop);
                        metapop.reset_rng(new_seed);
                        
                        new_seed = ea.rng().uniform_integer();
                        put<RNG_SEED>(new_seed, *knockout_loc2);
                        knockout_loc2->reset_rng(new_seed);
                        
                        if (mutations_off) {
                            put<TASK_MUTATION_PER_SITE_P>(0, metapop);
                            //put<MUTATION_PER_SITE_P>(0, metapop);
                            //put<GERM_MUTATION_PER_SITE_P>(0, metapop);
                        }
                        
                        typename EA::population_type init_mc;
                        init_mc.insert(init_mc.end(), knockout_loc2);
                        
//                        knockout_loc2->traits()._founder = i->traits().founder();
                        knockout_loc2->traits()._founder = ea.make_individual(*knockout_loc2);

                        std::swap(metapop.population(), init_mc);
                        
                        add_event<mt_gls_propagule>(metapop);
                        add_event<subpopulation_founder_event>(metapop);

                        int max_size = 32;
                        int max_update = 50000;
                        cur_update = 0;
                        
                        while ((metapop.size() < max_size) &&
                               (cur_update < max_update)){
                            metapop.update();
                            ++cur_update;
                        }
                        
                        if (track_details) {
                        std::string new_filename = "mt_gls_detail_uni_" + std::to_string(num_uni) + "_" + std::to_string(nr) + ".dat";
                        std::rename("mt_gls_detail.dat", new_filename.c_str());
                        }
                        
                        // get workload
                        float total_workload = 0;
                        int total_cells = 0;
                        int multicell_count = 0;
                        typedef typename EA::subpopulation_type::population_type subpop_type;
                                                
                        for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                            total_cells += j->size();
                            int cell_count = 0;
                            int org_id = get<ORGANISM_ID>(*j,0);
                            int org_parent_id = get<ORGANISM_PARENT_ID>(*j, 0);
                            

                            
                            for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                typename EA::subpopulation_type::individual_type& org=**m;
                                total_workload += get<WORKLOAD>(org, 0.0);
                                if (track_details) {
       
                                    df4.write(timepoint)
                                    .write(z)
                                    .write(q)
                                    .write(nr)
                                    .write(multicell_count)
                                    .write(org_id)
                                    .write(org_parent_id)
                                    .write(cell_count)
                                    .write(get<WORKLOAD>(org, 0.0))
                                    .write(get<GERM_STATUS>(org,true))
                                    .write(get<SAVED_RESOURCES>(org,0.0))
                                    .write(get<TASK_NOT>(org,0))
                                    .write(get<TASK_NAND>(org,0))
                                    .write(get<TASK_AND>(org,0))
                                    .write(get<TASK_ORNOT>(org,0))
                                    .write(get<TASK_OR>(org,0))
                                    .write(get<TASK_ANDNOT>(org,0))
                                    .write(get<TASK_NOR>(org,0))
                                    .write(get<TASK_XOR>(org,0))
                                    .write(get<TASK_EQUALS>(org,0))
                                    .write(cur_update)
                                    .write(get<GROUP_RESOURCE_UNITS>(*j,0.0))
                                    .write(total_workload);
                                    df4.write("\"");
                                    for(typename EA::subpopulation_type::genome_type::iterator k2=org.genome().begin(); k2!=org.genome().end(); ++k2) {
                                        df4.write(*k2)
                                        .write(" ");
                                    }
                                    df4.write("\"");
                                    df4.endl();

                                    
                                }
                                cell_count++;
                            }
                            multicell_count++;
                        }
                        
                        

                        
                        df.write(timepoint)
                        .write("uni")
                        .write(num_uni)
                        .write(nr)
                        .write(cur_update)
                        .write(total_workload)
                        .write(metapop.size())
                        .write(total_cells);
                        df.endl();
                        
                        
                    }
                    uni_count++;
                    
                }
                
            }
            
        }
        df2.write(timepoint)
        .write(num_uni)
        .write(num_uni_viable)
        .write(num_uni_inviable)
        .write(birth_up)
        .endl();
    }
    

    
    
        LIBEA_ANALYSIS_TOOL(lod_entrench) {
            
            datafile df("lod_entrench_all.dat");
            df.add_field("cost")
            .add_field("iteration")
            .add_field("update")
            .add_field("organism_size")
            .add_field("num_germ")
            .add_field("generation")
            .add_field("generation_diff")
            .add_field("workload")
            .add_field("workload_propagule_ineligible")
            ;
            
            
            datafile df2("lod_entrench_final.dat");
            df2.add_field("cost")
            .add_field("iteration")
            .add_field("update")
            .add_field("organism_size")
            .add_field("num_germ")
            .add_field("generation")
            .add_field("generation_diff")
            .add_field("workload")
            .add_field("workload_propagule_ineligible")
            .add_field("reverted")
            ;
            
            int num_rep = get<ANALYSIS_LOD_REPS>(ea,1);
            int mult = get<TISSUE_ACCRETION_MULT>(ea,0);
            int start_mult = mult - 6;
            int end_mult = mult + 6;
            if (start_mult < 0) {
                start_mult = 0;
            }
            int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);
            
            int meta_size = 1000;
            int entrench_not_found = true;
            std::set<int> checked_nums;
            
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i;
            if (timepoint == 1) {
                i = lod.end(); --i;
            } else {
                i=lod.begin(); i++;
                // find the first to transition
                for( ; i!=lod.end(); i++) {
                    if (i->size() > 2) {
                        break;
                    }
                }
            }
            
            while (entrench_not_found) {
                int revert_count = 0;
                checked_nums.insert(start_mult);
                for (int nr = 0; nr < num_rep; nr++) {
                    // should define checkpoint + analysis input
                    //ea is the thing loaded from the checkpoint; EA is its type
                    EA metapop; // a new EA
                    
                    metapop.initialize(ea.md());
                    put<TISSUE_ACCRETION_MULT>(start_mult, metapop);
                    
                    int new_seed = ea.rng().uniform_integer();
                    put<RNG_SEED>(new_seed, metapop);
                    metapop.reset_rng(new_seed);
                    
                    float start_gen = 0;
                    typename EA::population_type init_mc;
                    for (int j=0; j<meta_size; ++j){
                        typename EA::individual_ptr_type control_mc = metapop.make_individual(*i->traits().founder());
                        control_mc->initialize(metapop.md());
                        put<TISSUE_ACCRETION_MULT>(start_mult, *control_mc);
                        control_mc->reset_rng(metapop.rng().uniform_integer());
                        init_mc.insert(init_mc.end(),metapop.make_individual(*control_mc));
                        if (j ==0) {
                            start_gen = get<IND_GENERATION>(*control_mc);
                        }
                    }
                    
                    std::swap(metapop.population(), init_mc);
                    
                    add_event<mt_gls_propagule>(metapop);
                    
                    int max_update = 200000;
                    int cur_update = 0;
                    int exit = false;
                    
                    int exit_mean_size = 0;
                    
                    
                    while ((exit == false) &&
                           (cur_update < max_update)){
                        metapop.update();
                        ++cur_update;
                        
                        if ((cur_update % 100)==0) {
                            
                            float total_workload = 0;
                            float germ_workload = 0;
                            float organism_size = 0;
                            float num_germ = 0;
                            float gen = 0;
                            
                            typedef typename EA::subpopulation_type::population_type subpop_type;
                            
                            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                    typename EA::subpopulation_type::individual_type& org=**m;
                                    total_workload += get<WORKLOAD>(org, 0.0);
                                    if (get<GERM_STATUS>(org, 1)) {
                                        germ_workload += get<WORKLOAD>(org, 0.0);
                                        num_germ += 1;
                                    }
                                }
                                organism_size += j->population().size();
                                gen += get<IND_GENERATION>(*j);
                            }
                            
                            float mean_gen = gen/metapop.size();
                            float mean_gen_diff = mean_gen - start_gen;

                            float mean_size = organism_size/metapop.size();
                            if (metapop.current_update() > 5000){
                            if (mean_size < 2) {
                                exit_mean_size++;
                            } else {
                                exit_mean_size = 0;
                            }}
                            
                            df.write(start_mult)
                            .write(nr)
                            .write(metapop.current_update())
                            .write(organism_size/metapop.size())
                            .write(num_germ/metapop.size())
                            .write(mean_gen)
                            .write(mean_gen_diff)
                            .write(total_workload/organism_size)
                            .write(germ_workload/num_germ)
                            .endl();
                            
                            if ((exit_mean_size > 5) ||
                                (mean_gen_diff > 100))  {
                                int reverted = 0;
                                
                                if (exit_mean_size > 5)  {
                                    revert_count += 1;
                                    reverted = 1;
                                }
                                exit = true;
                                df2.write(start_mult)
                                .write(nr)
                                .write(metapop.current_update())
                                .write(organism_size/metapop.size())
                                .write(num_germ/metapop.size())
                                .write(mean_gen)
                                .write(mean_gen_diff)
                                .write(total_workload/organism_size)
                                .write(germ_workload/num_germ)
                                .write(reverted)
                                .endl();
                            }
                            if (cur_update == max_update){
                                df2.write(start_mult)
                                .write(nr)
                                .write(metapop.current_update())
                                .write(organism_size/metapop.size())
                                .write(num_germ/metapop.size())
                                .write(mean_gen)
                                .write(mean_gen_diff)
                                .write(total_workload/organism_size)
                                .write(germ_workload/num_germ)
                                .write("2")
                                .endl();
                            }
                        }
                    }
                    
                }
                if (start_mult <= end_mult) {
                    start_mult += 2;
                } else {
                    entrench_not_found = false;
                }

            }// end while
        }

    LIBEA_ANALYSIS_TOOL(lod_size) {
    
        datafile df("lod_size.dat");
        df.add_field("time")
        .add_field("organism_size")
        .add_field("num_germ")
        .add_field("workload")
        .add_field("germ_workload")
        ;
        
            
        line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);

        int i_count = 0;
        typename line_of_descent<EA>::iterator i;
        i=lod.begin(); i++;
        for( ; i!=lod.end(); i++) {

            typedef typename EA::subpopulation_type::population_type subpop_type;
            float total_workload = 0;
            float germ_workload = 0;
            int num_germ = 0;
            for(typename subpop_type::iterator m=i->population().begin(); m!=i->population().end(); ++m) {
                        typename EA::subpopulation_type::individual_type& org=**m;
                        total_workload += get<WORKLOAD>(org, 0.0);
                        if (get<GERM_STATUS>(org, 1)) {
                            germ_workload += get<WORKLOAD>(org, 0.0);
                            num_germ += 1;
                        }
            }

            df.write(i_count)
            .write(i->size())
            .write(num_germ)
            .write(total_workload)
            .write(germ_workload)
            .endl();
            
            i_count++;
        }
    }
    
   
    
    LIBEA_ANALYSIS_TOOL(lod_entrench_add) {
            
            datafile df("lod_entrench_all.dat");
            df.add_field("cost")
            .add_field("iteration")
            .add_field("update")
            .add_field("organism_size")
            .add_field("num_germ")
            .add_field("generation")
            .add_field("generation_diff")
            .add_field("workload")
            .add_field("workload_propagule_ineligible")
            ;
            
            
            datafile df2("lod_entrench_final.dat");
            df2.add_field("cost")
            .add_field("iteration")
            .add_field("update")
            .add_field("organism_size")
            .add_field("num_germ")
            .add_field("generation")
            .add_field("generation_diff")
            .add_field("workload")
            .add_field("workload_propagule_ineligible")
            .add_field("reverted")
            ;
            
            int num_rep = get<ANALYSIS_LOD_REPS>(ea,1);
            int add_ent = get<TISSUE_ACCRETION_ADD>(ea,0);
            int start_mult = add_ent;
            int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);
            
            int meta_size = 1000;
            int entrench_not_found = true;
            std::set<int> checked_nums;
            
            int i_count = 0;
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i; 
            if (timepoint == 1) {
                i = lod.end(); --i;
            } else {
                i=lod.begin(); i++;
                i_count++;
                // find the first to transition
                for( ; i!=lod.end(); i++) {
                    if (i->size() > 2) {
                        break;
                    }
                    i_count++;
                }
            }
            
            while (entrench_not_found) {
                int revert_count = 0;
                checked_nums.insert(start_mult);
                for (int nr = 0; nr < num_rep; nr++) {
                    // should define checkpoint + analysis input
                    //ea is the thing loaded from the checkpoint; EA is its type
                    EA metapop; // a new EA
                    
                    metapop.initialize(ea.md());
                    put<TISSUE_ACCRETION_ADD>(add_ent, metapop);
                    
                    int new_seed = ea.rng().uniform_integer();
                    put<RNG_SEED>(new_seed, metapop);
                    metapop.reset_rng(new_seed);
                    
                    float start_gen = 0;
                    typename EA::population_type init_mc;
                    for (int j=0; j<meta_size; ++j){
                        typename EA::individual_ptr_type control_mc = metapop.make_individual(*i->traits().founder());
                        control_mc->initialize(metapop.md());
                        put<TISSUE_ACCRETION_ADD>(start_mult, *control_mc);
                        control_mc->reset_rng(metapop.rng().uniform_integer());
                        init_mc.insert(init_mc.end(),metapop.make_individual(*control_mc));
                        if (j ==0) {
                            start_gen = get<IND_GENERATION>(*control_mc);
                        }
                    }
                    
                    std::swap(metapop.population(), init_mc);
                    
                    add_event<mt_gls_propagule>(metapop);
                    
                    int max_update = 200000;
                    int cur_update = 0;
                    int exit = false;
                    
                    int exit_mean_size = 0;
                    
                    
                    while ((exit == false) &&
                           (cur_update < max_update)){
                        metapop.update();
                        ++cur_update;
                        
                        if ((cur_update % 100)==0) {
                            
                            float total_workload = 0;
                            float germ_workload = 0;
                            float organism_size = 0;
                            float num_germ = 0;
                            float gen = 0;
                            
                            typedef typename EA::subpopulation_type::population_type subpop_type;
                            
                            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                    typename EA::subpopulation_type::individual_type& org=**m;
                                    total_workload += get<WORKLOAD>(org, 0.0);
                                    if (get<GERM_STATUS>(org, 1)) {
                                        germ_workload += get<WORKLOAD>(org, 0.0);
                                        num_germ += 1;
                                    }
                                }
                                organism_size += j->population().size();
                                gen += get<IND_GENERATION>(*j);
                            }
                            
                            float mean_gen = gen/metapop.size();
                            float mean_gen_diff = mean_gen - start_gen;

                            float mean_size = organism_size/metapop.size();
                            if (metapop.current_update() > 10000){
                            if (mean_size < 2) {
                                exit_mean_size++;
                            } else {
                                exit_mean_size = 0;
                            }}
                            
                            df.write(add_ent)
                            .write(nr)
                            .write(metapop.current_update())
                            .write(organism_size/metapop.size())
                            .write(num_germ/metapop.size())
                            .write(mean_gen)
                            .write(mean_gen_diff)
                            .write(total_workload/organism_size)
                            .write(germ_workload/num_germ)
                            .endl();
                            
                            if (mean_gen_diff > 100)  {
                                int reverted = 0;
                                
                                if (exit_mean_size > 5)  {
                                    revert_count += 1;
                                    reverted = 1;
                                }
                                exit = true;
                                df2.write(add_ent)
                                .write(nr)
                                .write(metapop.current_update())
                                .write(organism_size/metapop.size())
                                .write(num_germ/metapop.size())
                                .write(mean_gen)
                                .write(mean_gen_diff)
                                .write(total_workload/organism_size)
                                .write(germ_workload/num_germ)
                                .write(reverted)
                                .endl();
                            }
                            if (cur_update == max_update){
                                df2.write(add_ent)
                                .write(nr)
                                .write(metapop.current_update())
                                .write(organism_size/metapop.size())
                                .write(num_germ/metapop.size())
                                .write(mean_gen)
                                .write(mean_gen_diff)
                                .write(total_workload/organism_size)
                                .write(germ_workload/num_germ)
                                .write("2")
                                .endl();
                            }
                        }
                    }
                    
                }
                
                //add_ent += 100;
                add_ent *= 2;
                if (add_ent >= 4096) {
                    entrench_not_found = false;
                }
                
                
            }// end while
        }
    
    LIBEA_ANALYSIS_TOOL(lod_entrench_add_start_stop) {
        
        /* This one works backwards. Cost starts high and goes lower. Once we find a
         reversion point, then we stop. */
        
        datafile df("lod_entrench_all.dat");
        df.add_field("cost")
        .add_field("iteration")
        .add_field("update")
        .add_field("organism_size")
        .add_field("num_germ")
        .add_field("generation")
        .add_field("generation_diff")
        .add_field("workload")
        .add_field("workload_propagule_ineligible")
        ;
        
        
        datafile df2("lod_entrench_final.dat");
        df2.add_field("cost")
        .add_field("iteration")
        .add_field("update")
        .add_field("organism_size")
        .add_field("num_germ")
        .add_field("generation")
        .add_field("generation_diff")
        .add_field("workload")
        .add_field("workload_propagule_ineligible")
        .add_field("reverted")
        .add_field("lod")
        ;
        
        
        int num_rep = get<ANALYSIS_LOD_REPS>(ea,1);
        int lod_depth = 0;
        int next_lod = 0;
        
        
        int lod_start_analysis = get<LOD_START_ANALYSIS>(ea);
        int lod_end_analysis = get<LOD_END_ANALYSIS>(ea);
        
        line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
        typename line_of_descent<EA>::iterator i=lod.begin(); i++; i++;
        

        int add_ent = get<TISSUE_ACCRETION_ADD>(ea,0);
        int start_mult = add_ent;
        
        int meta_size = 1000;
        int entrench_not_found = true;
        std::set<int> checked_nums;
        
        
        for( ; i!=lod.end(); i++) {
            
                if (lod_depth != next_lod) {
                    lod_depth++;
                    continue;
                }
                
                next_lod += 100;
                
                if ((lod_depth < lod_start_analysis) ||
                    (lod_depth >= lod_end_analysis)) {
                    continue;
                }
            
            
        while (entrench_not_found) {
            int revert_count = 0;
            checked_nums.insert(start_mult);
            for (int nr = 0; nr < num_rep; nr++) {
                // should define checkpoint + analysis input
                //ea is the thing loaded from the checkpoint; EA is its type
                EA metapop; // a new EA
                
                metapop.initialize(ea.md());
                put<TISSUE_ACCRETION_ADD>(add_ent, metapop);
                
                int new_seed = ea.rng().uniform_integer();
                put<RNG_SEED>(new_seed, metapop);
                metapop.reset_rng(new_seed);
                
                float start_gen = 0;
                typename EA::population_type init_mc;
                for (int j=0; j<meta_size; ++j){
                    typename EA::individual_ptr_type control_mc = metapop.make_individual(*i->traits().founder());
                    control_mc->initialize(metapop.md());
                    put<TISSUE_ACCRETION_ADD>(start_mult, *control_mc);
                    control_mc->reset_rng(metapop.rng().uniform_integer());
                    init_mc.insert(init_mc.end(),metapop.make_individual(*control_mc));
                    if (j ==0) {
                        start_gen = get<IND_GENERATION>(*control_mc);
                    }
                }
                
                std::swap(metapop.population(), init_mc);
                
                add_event<mt_gls_propagule>(metapop);
                
                int max_update = 200000;
                int cur_update = 0;
                int exit = false;
                
                int exit_mean_size = 0;
                
                
                while ((exit == false) &&
                       (cur_update < max_update)){
                    metapop.update();
                    ++cur_update;
                    
                    if ((cur_update % 100)==0) {
                        
                        float total_workload = 0;
                        float germ_workload = 0;
                        float organism_size = 0;
                        float num_germ = 0;
                        float gen = 0;
                        
                        typedef typename EA::subpopulation_type::population_type subpop_type;
                        
                        for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                            for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                                typename EA::subpopulation_type::individual_type& org=**m;
                                total_workload += get<WORKLOAD>(org, 0.0);
                                if (get<GERM_STATUS>(org, 1)) {
                                    germ_workload += get<WORKLOAD>(org, 0.0);
                                    num_germ += 1;
                                }
                            }
                            organism_size += j->population().size();
                            gen += get<IND_GENERATION>(*j);
                        }
                        
                        float mean_gen = gen/metapop.size();
                        float mean_gen_diff = mean_gen - start_gen;

                        float mean_size = organism_size/metapop.size();
                        if (metapop.current_update() > 10000){
                        if (mean_size < 2) {
                            exit_mean_size++;
                        } else {
                            exit_mean_size = 0;
                        }}
                        
                        df.write(add_ent)
                        .write(nr)
                        .write(metapop.current_update())
                        .write(organism_size/metapop.size())
                        .write(num_germ/metapop.size())
                        .write(mean_gen)
                        .write(mean_gen_diff)
                        .write(total_workload/organism_size)
                        .write(germ_workload/num_germ)
                        .endl();
                        
                        if (mean_gen_diff > 100)  {
                            int reverted = 0;
                            
                            if (exit_mean_size > 5)  {
                                revert_count += 1;
                                reverted = 1;
                            }
                            exit = true;
                            df2.write(add_ent)
                            .write(nr)
                            .write(metapop.current_update())
                            .write(organism_size/metapop.size())
                            .write(num_germ/metapop.size())
                            .write(mean_gen)
                            .write(mean_gen_diff)
                            .write(total_workload/organism_size)
                            .write(germ_workload/num_germ)
                            .write(reverted)
                            .write(lod_depth)
                            .endl();
                        }
                        if (cur_update == max_update){
                            df2.write(add_ent)
                            .write(nr)
                            .write(metapop.current_update())
                            .write(organism_size/metapop.size())
                            .write(num_germ/metapop.size())
                            .write(mean_gen)
                            .write(mean_gen_diff)
                            .write(total_workload/organism_size)
                            .write(germ_workload/num_germ)
                            .write("2")
                            .write(lod_depth)
                            .endl();
                        }
                    }
                }
                
            }
            
            add_ent /= 2;
            if ((add_ent >= 2048) || (add_ent <= 1)) {
                entrench_not_found = false;
            }
            
            /*if (revert_count == 0) {
                entrench_not_found = false;
            }*/
            
            
        }// end while
        }// end for
    }
    
//    LIBEA_ANALYSIS_TOOL(lod_entrench_mc_cost) {
//
//            datafile df("lod_entrench_all.dat");
//            df.add_field("cost")
//            .add_field("iteration")
//            .add_field("update")
//            .add_field("organism_size")
//            .add_field("num_germ")
//            .add_field("generation")
//            .add_field("generation_diff")
//            .add_field("workload")
//            .add_field("workload_propagule_ineligible")
//            ;
//
//
//            datafile df2("lod_entrench_final.dat");
//            df2.add_field("cost")
//            .add_field("iteration")
//            .add_field("update")
//            .add_field("organism_size")
//            .add_field("num_germ")
//            .add_field("generation")
//            .add_field("generation_diff")
//            .add_field("workload")
//            .add_field("workload_propagule_ineligible")
//            .add_field("reverted")
//            ;
//
//            int num_rep = get<ANALYSIS_LOD_REPS>(ea,1);
//            int start_rep = get<ANALYSIS_LOD_START_REP>(ea,0);
//
//            int mc_res = get<GROUP_REP_THRESHOLD>(ea);
//            int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);
//
//            int meta_size = 1000;
//            int entrench_not_found = true;
//            std::set<int> checked_nums;
//
//
//            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
//            typename line_of_descent<EA>::iterator i;
//            if (timepoint == 1) {
//                i = lod.end(); --i;
//            } else {
//                i=lod.begin(); i++;
//                // find the first to transition
//                for( ; i!=lod.end(); i++) {
//                    if (i->size() > 2) {
//                        break;
//                    }
//                }
//            }
//
//            while (entrench_not_found) {
//                int revert_count = 0;
//                for (int nr = start_rep; nr < start_rep+num_rep; nr++) {
//                    // should define checkpoint + analysis input
//                    //ea is the thing loaded from the checkpoint; EA is its type
//                    EA metapop; // a new EA
//
//                    metapop.initialize(ea.md());
//                    put<GROUP_REP_THRESHOLD>(mc_res, metapop);
//
//                    int new_seed = ea.rng().uniform_integer();
//                    put<RNG_SEED>(new_seed, metapop);
//                    metapop.reset_rng(new_seed);
//
//                    float start_gen = 0;
//                    typename EA::population_type init_mc;
//                    for (int j=0; j<meta_size; ++j){
//                        typename EA::individual_ptr_type control_mc = metapop.make_individual(*i->traits().founder());
//                        control_mc->initialize(metapop.md());
//                        put<GROUP_REP_THRESHOLD>(mc_res, *control_mc);
//                        control_mc->reset_rng(metapop.rng().uniform_integer());
//                        init_mc.insert(init_mc.end(),metapop.make_individual(*control_mc));
//                        if (j ==0) {
//                            start_gen = get<IND_GENERATION>(*control_mc);
//                        }
//                    }
//
//                    std::swap(metapop.population(), init_mc);
//
//                    add_event<mt_gls_propagule>(metapop);
//
//                    int max_update = 200000;
//                    int cur_update = 0;
//                    int exit = false;
//
//                    int exit_mean_size = 0;
//
//
//                    while ((exit == false) &&
//                           (cur_update < max_update)){
//                        metapop.update();
//                        ++cur_update;
//
//                        if ((cur_update % 100)==0) {
//
//                            float total_workload = 0;
//                            float germ_workload = 0;
//                            float organism_size = 0;
//                            float num_germ = 0;
//                            float gen = 0;
//
//                            typedef typename EA::subpopulation_type::population_type subpop_type;
//
//                            for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
//                                for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
//                                    typename EA::subpopulation_type::individual_type& org=**m;
//                                    total_workload += get<WORKLOAD>(org, 0.0);
//                                    if (get<GERM_STATUS>(org, 1)) {
//                                        germ_workload += get<WORKLOAD>(org, 0.0);
//                                        num_germ += 1;
//                                    }
//                                }
//                                organism_size += j->population().size();
//                                gen += get<IND_GENERATION>(*j);
//                            }
//
//                            float mean_gen = gen/metapop.size();
//                            float mean_gen_diff = mean_gen - start_gen;
//
//                            float mean_size = organism_size/metapop.size();
//                            if (metapop.current_update() > 1000){
//                                if (mean_size < 2) {
//                                    exit_mean_size++;
//                                } else  {
//                                    exit_mean_size = 0;
//                            }}
//
//
//
//
//                            df.write(mc_res)
//                            .write(nr)
//                            .write(metapop.current_update())
//                            .write(organism_size/metapop.size())
//                            .write(num_germ/metapop.size())
//                            .write(mean_gen)
//                            .write(mean_gen_diff)
//                            .write(total_workload/organism_size)
//                            .write(germ_workload/num_germ)
//                            .endl();
//
//                            if ((exit_mean_size >= 5) ||
//                                (mean_gen_diff > 100))  {
//                                int reverted = 0;
//
//                                if ((exit_mean_size >= 5) &&
//                                    (mean_gen_diff > 10)){
//                                    revert_count += 1;
//                                    reverted = 1;
//                                }
//                                exit = true;
//                                df2.write(mc_res)
//                                .write(nr)
//                                .write(metapop.current_update())
//                                .write(organism_size/metapop.size())
//                                .write(num_germ/metapop.size())
//                                .write(mean_gen)
//                                .write(mean_gen_diff)
//                                .write(total_workload/organism_size)
//                                .write(germ_workload/num_germ)
//                                .write(reverted)
//                                .endl();
//                            }
//                            if (cur_update == max_update){
//                                df2.write(mc_res)
//                                .write(nr)
//                                .write(metapop.current_update())
//                                .write(organism_size/metapop.size())
//                                .write(num_germ/metapop.size())
//                                .write(mean_gen)
//                                .write(mean_gen_diff)
//                                .write(total_workload/organism_size)
//                                .write(germ_workload/num_germ)
//                                .write("2")
//                                .endl();
//                            }
//                        }
//                    }
//
//                }
//
//                mc_res -= 25;
//                if (mc_res < 0) {
//                    entrench_not_found = false;
//                }
//
//
//            }// end while
//        }
    
    
    
LIBEA_ANALYSIS_TOOL(lod_dol) {
    
    datafile df("lod_dol.dat");

    df.add_field("timepoint")
    .add_field("tech_rep")
    .add_field("num_propagules")
    .add_field("cur_update")
    .add_field("time_since_last_rep")
    .add_field("multicell_size")
    .add_field("germ_count")
    .add_field("total_workload")
    .add_field("germ_workload")
    .add_field("soma_workload");
    

    int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);

    // get right lod member
    line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
    typename line_of_descent<EA>::iterator i;
    int lod_step = 0;

    if (timepoint == 1) {
        i = lod.end(); --i;
    } else if (timepoint == 0) {
        i=lod.begin(); i++;
        lod_step++;
        // find the first to transition
        for( ; i!=lod.end(); i++) {
            lod_step++;
            if ((timepoint != 0) && (timepoint == lod_step)) {
                break;
            }
            if ((timepoint == 0) && (i->size() > 2)) {
                break;
            }
        }
    } else {
        i=lod.begin(); i++;
        lod_step++;
        for( ; i!=lod.end(); i++) {
            lod_step++;
            if (timepoint == lod_step) {
                break;
            }
        }
    }
         


    int num_rep = 100;
    for (int nr = 0; nr < num_rep; nr++) {

        typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
        control_ea->initialize(ea.md());
        control_ea->reset_rng(ea.rng().uniform_integer());

        // run until X replication events have 'happened'. Really, but happened, here we just add the resources back in.
        int cur_update = 0;
        int time_since_last_rep = 0;
        int update_max = 2000;
        int num_replications = 0;
        //
        while ((cur_update < update_max) && (num_replications < 5)) {
            control_ea->update();
            ++cur_update;
            ++time_since_last_rep;


            if (get<DIVIDE_REMOTE>(*control_ea,0)) {
                num_replications++;

                int total_workload = 0;
                int germ_workload = 0;
                int germ_count = 0;

                typedef typename EA::subpopulation_type::population_type subpop_type;

                for(typename subpop_type::iterator m=control_ea->population().begin(); m!=control_ea->population().end(); ++m) {
                    typename EA::subpopulation_type::individual_type& org=**m;
                    total_workload += get<WORKLOAD>(org, 0.0);
                    if (get<GERM_STATUS>(org, true)) {
                        ++germ_count;
                        germ_workload += get<WORKLOAD>(org, 0.0);
                    }
                }

                df.write(timepoint)
                .write(nr)
                .write(num_replications)
                .write(cur_update)
                .write(time_since_last_rep)
                .write(control_ea->size())
                .write(germ_count)
                .write(total_workload)
                .write(germ_workload)
                .write(total_workload - germ_workload);
                df.endl();


                control_ea->resources().reset();
                put<GROUP_RESOURCE_UNITS>(0,*control_ea);
                put<DIVIDE_REMOTE>(0,*control_ea);
                time_since_last_rep = 0;

            }
        }
            
    }// end for
}//end lod
    
  
    LIBEA_ANALYSIS_TOOL(lod_task_switching_dol) {
        
        datafile df("lod_dol.dat");

        df.add_field("timepoint")
        .add_field("tech_rep")
        .add_field("num_propagules")
        .add_field("cur_update")
        .add_field("time_since_last_rep")
        .add_field("multicell_size")
        .add_field("germ_count")
        .add_field("num_task_switches")
        .add_field("mean_task_switches")
        .add_field("not")
        .add_field("nand")
        .add_field("and")
        .add_field("ornot")
        .add_field("or")
        .add_field("andnot")
        .add_field("nor")
        .add_field("xor")
        .add_field("equals")
        .add_field("total_tasks")
        .add_field("dol_metric")
        .add_field("shannon_mutual_info")
        .add_field("shannon_mutual_info_norm");
        

        int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);

        // get right lod member
        line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
        typename line_of_descent<EA>::iterator i;
        int lod_step = 0;

        if (timepoint == 1) {
            i = lod.end(); --i;
        } else if (timepoint == 0) {
            i=lod.begin(); i++;
            lod_step++;
            // find the first to transition
            for( ; i!=lod.end(); i++) {
                lod_step++;
                if ((timepoint != 0) && (timepoint == lod_step)) {
                    break;
                }
                if ((timepoint == 0) && (i->size() > 2)) {
                    break;
                }
            }
        } else {
            i=lod.begin(); i++;
            lod_step++;
            for( ; i!=lod.end(); i++) {
                lod_step++;
                if (timepoint == lod_step) {
                    break;
                }
            }
        }
             


        int num_rep = 100;
        for (int nr = 0; nr < num_rep; nr++) {
            std::vector<std::string> tps;

            typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
            control_ea->initialize(ea.md());
            control_ea->reset_rng(ea.rng().uniform_integer());
            //add_event<task_profile_tracking>(*control_ea);
            //add_event<task_profile_birth_event>(*control_ea);

            // run until X replication events have 'happened'. Really, but happened, here we just add the resources back in.
            int cur_update = 0;
            int time_since_last_rep = 0;
            int update_max = 2000;
            int num_replications = 0;
            float ts = 0;
            //
            int t_not = 0;
            int t_nand = 0;
            int t_and = 0;
            int t_ornot = 0;
            int t_or = 0;
            int t_andnot = 0;
            int t_nor = 0;
            int t_xor = 0;
            int t_equals = 0;
            
            while ((cur_update < update_max) && (num_replications < 5)) {
                control_ea->update();
                ++cur_update;
                ++time_since_last_rep;


                if (get<DIVIDE_REMOTE>(*control_ea,0)) {
                    num_replications++;

                    int germ_count = 0;

                    typedef typename EA::subpopulation_type::population_type subpop_type;
                    std::vector< std::vector<double> > pij;
                    std::vector<double> pj (9);
                    double pop_count = 0;
                    double active_pop = 0;

                    for(typename subpop_type::iterator m=control_ea->population().begin(); m!=control_ea->population().end(); ++m) {
                        typename EA::subpopulation_type::individual_type& org=**m;


                        //tps.push_back(get<TASK_PROFILE>(org,""));
                        ts += get<NUM_SWITCHES>(org, 0);
                        
                        t_not += get<TASK_NOT>(org, 0.0);
                        t_nand += get<TASK_NAND>(org, 0.0);
                        t_and += get<TASK_AND>(org, 0.0);
                        t_ornot += get<TASK_ORNOT>(org, 0.0);
                        t_or += get<TASK_OR>(org, 0.0);
                        t_andnot += get<TASK_ANDNOT>(org, 0.0);
                        t_nor += get<TASK_NOR>(org, 0.0);
                        t_xor += get<TASK_XOR>(org, 0.0);
                        t_equals += get<TASK_EQUALS>(org, 0.0);

                        if (get<GERM_STATUS>(org, true)) {
                            ++germ_count;
                            
                        }
                    
                        ++pop_count;
                        std::vector<double> porg (9);
                        porg[0] = get<TASK_NOT>(org,0.0);
                        porg[1] = get<TASK_NAND>(org,0.0);
                        porg[2] = get<TASK_AND>(org,0.0);
                        porg[3] = get<TASK_ORNOT>(org,0.0);
                        porg[4] = get<TASK_OR>(org,0.0);
                        porg[5] = get<TASK_ANDNOT>(org,0.0);
                        porg[6] = get<TASK_NOR>(org,0.0);
                        porg[7] = get<TASK_XOR>(org,0.0);
                        porg[8] = get<TASK_EQUALS>(org,0.0);
                        
                        double total_num_tasks = std::accumulate(porg.begin(), porg.end(), 0);
                        
                        // Normalize the tasks and add to matrix
                        if(total_num_tasks > 0) {
                            for (unsigned int k=0; k<porg.size(); ++k) {
                                porg[k] /= total_num_tasks;
                            }
                            ++active_pop;
                            pij.push_back(porg);
                        }
                    }
                        
                    double shannon_sum = 0.0;
                    double shannon_norm = 0.0;
                    if (active_pop > 1) {
                        // figure out pj
                        for (unsigned int k=0; k<pj.size(); ++k) {
                            for (int m=0; m<active_pop; ++m) {
                                pj[k] += pij[m][k];
                            }
                            pj[k] /= active_pop;
                        }
                        
                        // compute shannon mutual information based on matrix...
                        double shannon_change = 0.0;
                        double t_pij = 0.0;
                        double t_pi = 1.0/active_pop;
                        double t_pj = 0;
                        double pij_sum = 0.0;
                        // calculate shannon mutual information
                        for (unsigned int i=0; i<active_pop; i++) {
                            for (int j=0; j<pj.size(); j++) {
                                t_pij = pij[i][j]/active_pop;
                                t_pj = pj[j];
                                pij_sum += t_pij;
                                if (t_pi && t_pj && t_pij) {
                                    shannon_change= (t_pij * log(t_pij / (t_pi * t_pj)));
                                    shannon_sum += shannon_change;
                                }
                            }
                        }
                    }
                    if ((shannon_sum > 0 ) &&(active_pop > 0)) {
                        shannon_norm = shannon_sum / log((double)active_pop);
                    }
                        
//                        shannon_sum_all += shannon_sum;
//                        shannon_norm_all += shannon_norm;
//                        active_pop_all += active_pop;
//                        pop_count_all += pop_count;
//
//
//                    }
//                    _df.write(shannon_sum_all / num_multis)
//                    .write(shannon_norm_all / num_multis)
//                    .write(active_pop_all / num_multis)
//                    .write(pop_count_all / num_multis);
//                    _df.endl();
            
            
                    
                   
                    
                    float mean_task_switches = ts/control_ea->size();
                    float total_tasks = t_not + t_nand + t_and + t_ornot + t_nor + t_xor + t_equals + t_andnot + t_or;

                    df.write(timepoint)
                    .write(nr)
                    .write(num_replications)
                    .write(cur_update)
                    .write(time_since_last_rep)
                    .write(control_ea->size())
                    .write(germ_count)
                    .write(ts)
                    .write(mean_task_switches)
                    .write(t_not)
                    .write(t_nand)
                    .write(t_and)
                    .write(t_ornot)
                    .write(t_or)
                    .write(t_andnot)
                    .write(t_nor)
                    .write(t_xor)
                    .write(t_equals)
                    .write(total_tasks)
                    .write(total_tasks/mean_task_switches)
                    .write(shannon_sum)
                    .write(shannon_norm);
                    df.endl();


                    control_ea->resources().reset();
                    put<GROUP_RESOURCE_UNITS>(0,*control_ea);
                    put<DIVIDE_REMOTE>(0,*control_ea);
                    time_since_last_rep = 0;

                }
            }
                
        }// end for
 
    }//end lod
    
    
    LIBEA_ANALYSIS_TOOL(lod_report_gs) {
        
        line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
        
        typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
        
        datafile df("lod_report_gs.dat");
        df.add_field("lod_depth")
        .add_field("update")
        .add_field("fit")
        .add_field("size")
        .add_field("num_germ")
        .add_field("num_soma")
        .add_field("germ_workload")
        .add_field("germ_workload_var")
        .add_field("soma_workload")
        .add_field("soma_workload_var")
        ;
        
        
        datafile df2("lod_fitness.dat");
        df2.add_field("lod")
        .add_field("iteration")
        .add_field("time_to_fill")
        .add_field("workload")
        .add_field("num_org")
        ;
        
        
        int num_rep = get<ANALYSIS_LOD_REPS>(ea,1);

        
        int lod_depth = 0;
        // skip def ancestor (that's what the +1 does)
        for( ; i!=lod.end(); ++i) {
            
            df.write(lod_depth);
            df.write(get<IND_BIRTH_UPDATE>(*i->traits().founder(),0));

            // **i is the EA, AS OF THE TIME THAT IT DIED!
            
            // To replay, need to create new eas for each knockout exper.
            // setup the population (really, an ea):
            typename EA::individual_ptr_type control_ea = ea.make_individual(*i->traits().founder());
            
            // replay! till the group amasses the right amount of resources
            // or exceeds its window...
            int cur_update = 0;
            int update_max = 2000;
            
            // and run till the group amasses the right amount of resources
            while ((get<DIVIDE_REMOTE>(*control_ea,0) == 0) &&
                   (cur_update < update_max)){
                control_ea->update();
                ++cur_update;
            }
            
            df.write(cur_update);
            df.write(control_ea->population().size());
            
            int germ_count = 0;
            int soma_count = 0;
            accumulator_set<double, stats<tag::mean, tag::variance> > germ_workload_acc;
            accumulator_set<double, stats<tag::mean, tag::variance> > soma_workload_acc;
            
            for(typename EA::subpopulation_type::population_type::iterator j=control_ea->population().begin(); j!=control_ea->population().end(); ++j) {
                
                typename EA::subpopulation_type::individual_type& org=**j;
                if (get<GERM_STATUS>(org, true)) {
                    ++germ_count;
                    germ_workload_acc(get<WORKLOAD>(org, 0.0));
                } else {
                    soma_workload_acc(get<WORKLOAD>(org, 0.0));
                    ++soma_count;
                }
            }
            
            
            if (germ_count) {
                df.write(germ_count)
                .write(soma_count)
                .write(mean(germ_workload_acc))
                .write(variance(germ_workload_acc))
                ;
            } else {
                df.write(0)
                .write(0)
                .write(0)
                .write(0);
            }
            if (soma_count) {
            df.write(mean(soma_workload_acc))
            .write(variance(soma_workload_acc));
            } else {
                df.write(0)
                .write(0);
            }
            
            for (int nr = 0; nr < num_rep; nr++) {
                
                
                // should define checkpoint + analysis input
                //ea is the thing loaded from the checkpoint; EA is its type
                EA metapop; // a new EA
                //typedef metadata md_type;
                
                typename EA::md_type md(ea.md());
                // override md settings (pop size, geo, etc)
                
                
                metapop.initialize(md);
                put<METAPOPULATION_SIZE>(32, metapop);
                put<RUN_UPDATES>(10000, metapop);
                put<RNG_SEED>(nr, metapop);
                
                if (nr != 0) {
                    metapop.reset_rng(nr);
                }
                
                typename EA::individual_ptr_type control_mc = ea.make_individual(*i->traits().founder());
                put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *control_mc);
                
                typename EA::population_type init_mc;
                init_mc.insert(init_mc.end(),control_mc);
                
                std::swap(metapop.population(), init_mc);
                
                add_event<mt_gls_propagule>(metapop);
                
                int max_size = 32;
                int max_update = 50000;
                int cur_update = 0;
                
                while ((metapop.size() < max_size) &&
                       (cur_update < max_update)){
                    metapop.update();
                    ++cur_update;
                }
                
                // get workload
                float total_workload = 0;
                typedef typename EA::subpopulation_type::population_type subpop_type;
                
                for(typename EA::iterator j=metapop.begin(); j!=metapop.end(); ++j) {
                    for(typename subpop_type::iterator m=j->population().begin(); m!=j->population().end(); ++m) {
                        typename EA::subpopulation_type::individual_type& org=**m;
                        total_workload += get<WORKLOAD>(org, 0.0);
                    }
                }
                
                df2.write(lod_depth)
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df2.endl();
            }

            
            
            df.endl();
            
            ++lod_depth;
        }
        
    }
        
    }
}


#endif
