
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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>


LIBEA_MD_DECL(LOD_START_ANALYSIS, "ea.mt.lod_start_analysis", int);
LIBEA_MD_DECL(LOD_END_ANALYSIS, "ea.mt.lod_end_analysis", int);
LIBEA_MD_DECL(ANALYSIS_LOD_REPS, "ea.mt.lod_analysis_reps", int);
LIBEA_MD_DECL(ANALYSIS_LOD_START_COST, "ea.mt.lod_start_cost", int);
LIBEA_MD_DECL(ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE, "ea.mt.lod_timepoint_to_analyze", int);
LIBEA_MD_DECL(ANALYSIS_LOD_START_REP, "ea.mt.lod_analysis_start_rep", int);




namespace ealib {
    namespace analysis {
        
        
        
        
        LIBEA_ANALYSIS_TOOL(lod_fitness) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;

            int num_rep = 100;

            //lifecycle::after_initialization(metapop);
            //lifecycle::gather_events(metapop);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.end(); --i;
            
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
            
            
                df.write("final")
                .write("mc")
                .write("0")
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df.endl();
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


            for (int z =0; z < 100; z++) {
                for (int q = 0; q < control_ea->isa().size(); q++) {

                    typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                    typename EA::individual_ptr_type knockout_loc2 = ea.make_individual(*i->traits().founder());
                    
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc2);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc2);
                    
                    knockout_loc->population()[0]->genome()[z] = q;
                    knockout_loc2->population()[0]->genome()[z] = q;
                    
                    put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    
                    
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                           (cur_update < update_max)){
                        knockout_loc->update();
                        ++cur_update;
                    }
                    


                    
                    if (knockout_loc->population().size() < 2) {
                        num_uni++;
                        
                        if (cur_update == update_max) {
                            num_uni_inviable++;
                            continue;
                        }
                        
                        num_uni_viable++;
                        

                        for (int nr = 0; nr < num_rep; nr++) {
                            // should define checkpoint + analysis input
                            //ea is the thing loaded from the checkpoint; EA is its type
                            EA metapop; // a new EA
                            //typedef metadata md_type;
                            typename EA::md_type md(ea.md());
                            // override md settings (pop size, geo, etc)
                            
                            
                            metapop.initialize(md);
                            put<METAPOPULATION_SIZE>(32, metapop);
                            put<RUN_UPDATES>(50000, metapop);
                            put<RNG_SEED>(nr, metapop);
                            if (nr != 0) {
                                metapop.reset_rng(nr);
                            }
                            
                            typename EA::population_type init_mc;
                            init_mc.insert(init_mc.end(),knockout_loc2);
                            
                            std::swap(metapop.population(), init_mc);
                            
                            add_event<mt_gls_propagule>(metapop);
                            
                            int max_size = 32;
                            int max_update = 50000;
                            cur_update = 0;
                            
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
                            
                            df.write("final")
                            .write("uni")
                            .write(uni_count)
                            .write(nr)
                            .write(cur_update)
                            .write(total_workload)
                            .write(metapop.size());
                            df.endl();
                            
                            
                        }
                        uni_count++;

                    }
                    
                }
                
            }
            df2.write("final")
            .write(num_uni)
            .write(num_uni_viable)
            .write(num_uni_inviable)
            .write(birth_up)
            .endl();

            
        }
        
        LIBEA_ANALYSIS_TOOL(lod_fitness_no_mutations) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("update")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;
            
            int num_rep = 100;
            
            //lifecycle::after_initialization(metapop);
            //lifecycle::gather_events(metapop);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.end(); --i;
            
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
                put<TASK_MUTATION_PER_SITE_P>(0, metapop);
                put<MUTATION_PER_SITE_P>(0, metapop);
                put<GERM_MUTATION_PER_SITE_P>(0, metapop);

                
                if (nr != 0) {
                    metapop.reset_rng(nr);
                }
                
                typename EA::individual_ptr_type control_mc = ea.make_individual(*i->traits().founder());
                int birth_up = get<IND_BIRTH_UPDATE>(*i->traits().founder(),0);

                put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *control_mc);
                put<TASK_MUTATION_PER_SITE_P>(0, *control_mc);
                put<MUTATION_PER_SITE_P>(0, *control_mc);

                
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
                
                
                df.write("final")
                .write("mc")
                .write("0")
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df.endl();
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
            
            
            for (int z =0; z < 100; z++) {
                for (int q = 0; q < control_ea->isa().size(); q++) {
                    
                    typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                    typename EA::individual_ptr_type knockout_loc2 = ea.make_individual(*i->traits().founder());
                    
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc2);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc2);
                    
                    knockout_loc->population()[0]->genome()[z] = q;
                    knockout_loc2->population()[0]->genome()[z] = q;
                    
                    put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    
                    put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc2);
                    put<MUTATION_PER_SITE_P>(0, *knockout_loc2);
                    put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc2);

                    
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                           (cur_update < update_max)){
                        knockout_loc->update();
                        ++cur_update;
                    }
                    
                    
                    
                    
                    if (knockout_loc->population().size() < 2) {
                        num_uni++;
                        
                        if (cur_update == update_max) {
                            num_uni_inviable++;
                            continue;
                        }
                        
                        num_uni_viable++;
                        
                        
                        for (int nr = 0; nr < num_rep; nr++) {
                            // should define checkpoint + analysis input
                            //ea is the thing loaded from the checkpoint; EA is its type
                            EA metapop; // a new EA
                            //typedef metadata md_type;
                            typename EA::md_type md(ea.md());
                            // override md settings (pop size, geo, etc)
                            
                            
                            metapop.initialize(md);
                            put<METAPOPULATION_SIZE>(32, metapop);
                            put<RUN_UPDATES>(50000, metapop);
                            put<RNG_SEED>(nr, metapop);
                            put<TASK_MUTATION_PER_SITE_P>(0, metapop);
                            put<MUTATION_PER_SITE_P>(0, metapop);
                            put<GERM_MUTATION_PER_SITE_P>(0, metapop);

                            if (nr != 0) {
                                metapop.reset_rng(nr);
                            }
                            
                            typename EA::population_type init_mc;
                            init_mc.insert(init_mc.end(),knockout_loc2);
                            
                            std::swap(metapop.population(), init_mc);
                            
                            add_event<mt_gls_propagule>(metapop);
                            
                            int max_size = 32;
                            int max_update = 50000;
                            cur_update = 0;
                            
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
                            
                            df.write("final")
                            .write("uni")
                            .write(uni_count)
                            .write(nr)
                            .write(cur_update)
                            .write(total_workload)
                            .write(metapop.size());
                            df.endl();
                            
                            
                        }
                        uni_count++;
                        
                    }
                    
                }
                
            }
            df2.write("final")
            .write(num_uni)
            .write(num_uni_viable)
            .write(num_uni_inviable)
            .write(birth_up)
            .endl();
            
            
        }

        
        LIBEA_ANALYSIS_TOOL(lod_fitness_at_trans) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;
            
            int num_rep = 100;
            
            //lifecycle::after_initialization(metapop);
            //lifecycle::gather_events(metapop);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.begin(); i++;
            
            // find the first to transition
            for( ; i!=lod.end(); i++) {
                if (i->size() > 2) {
                    break;
                }
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
                
                
                df.write("trans")
                .write("mc")
                .write("0")
                .write(nr)
                .write(cur_update)
                .write(total_workload)
                .write(metapop.size());
                df.endl();
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
            
            
            for (int z =0; z < 100; z++) {
                for (int q = 0; q < control_ea->isa().size(); q++) {
                    typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                    typename EA::individual_ptr_type knockout_loc2 = ea.make_individual(*i->traits().founder());
                    
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                    put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc2);
                    put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc2);
                    
                    knockout_loc->population()[0]->genome()[z] = q;
                    knockout_loc2->population()[0]->genome()[z] = q;

                    put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<MUTATION_PER_SITE_P>(0, *knockout_loc);
                    put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc);
                    
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                           (cur_update < update_max)){
                        knockout_loc->update();
                        ++cur_update;
                    }

                    
                    
                    if (knockout_loc->population().size() < 2) {
                        num_uni++;
                        
                        if (cur_update == update_max) {
                            num_uni_inviable++;
                            continue;
                        }

                        
                        num_uni_viable++;
                        
                        
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
                            
                            typename EA::population_type init_mc;
                            init_mc.insert(init_mc.end(),knockout_loc2);
                            
                            std::swap(metapop.population(), init_mc);
                            
                            add_event<mt_gls_propagule>(metapop);
                            
                            int max_size = 32;
                            int max_update = 50000;
                            cur_update = 0;
                            
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
                            
                            df.write("trans")
                            .write("uni")
                            .write(num_uni)
                            .write(nr)
                            .write(cur_update)
                            .write(total_workload)
                            .write(metapop.size());
                            df.endl();
                            
                            
                        }
                        uni_count++;
                        
                    }
                    
                }
                
            }
            df2.write("trans")
            .write(num_uni)
            .write(num_uni_viable)
            .write(num_uni_inviable)
            .write(birth_up)
            .endl();
        }
        
        
        LIBEA_ANALYSIS_TOOL(lod_fitness_start_stop) {
            
            datafile df("lod_fitness.dat");
            df.add_field("timepoint")
            .add_field("mc_or_uni")
            .add_field("count")
            .add_field("iteration")
            .add_field("time_to_fill")
            .add_field("workload")
            .add_field("num_org")
            ;
            
            
            datafile df2("lod_fit_summary.dat");
            df2.add_field("timepoint")
            .add_field("num_unicell_revertants")
            .add_field("num_viable_unicells")
            .add_field("num_inviable_unicells")
            .add_field("update")
            ;
            
            int num_rep = 100;
            
            int lod_depth = 0;
            int next_lod = 0;
            
            
            int lod_start_analysis = get<LOD_START_ANALYSIS>(ea);
            int lod_end_analysis = get<LOD_END_ANALYSIS>(ea);
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            typename line_of_descent<EA>::iterator i=lod.begin(); i++; i++;
            
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
                        
                        
                        df.write(lod_depth)
                        .write("mc")
                        .write("0")
                        .write(nr)
                        .write(cur_update)
                        .write(total_workload)
                        .write(metapop.size());
                        df.endl();
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
                    
                    
                    for (int z =0; z < 100; z++) {
                        for (int q = 0; q < control_ea->isa().size(); q++) {
                            typename EA::individual_ptr_type knockout_loc = ea.make_individual(*i->traits().founder());
                            typename EA::individual_ptr_type knockout_loc2 = ea.make_individual(*i->traits().founder());
                            
                            put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc);
                            put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc);
                            put<IND_REP_THRESHOLD>(get<IND_REP_THRESHOLD>(ea,0), *knockout_loc2);
                            put<COST_START_UPDATE>(get<COST_START_UPDATE>(ea,0), *knockout_loc2);
                            
                            knockout_loc->population()[0]->genome()[z] = q;
                            knockout_loc2->population()[0]->genome()[z] = q;
                            
                            put<TASK_MUTATION_PER_SITE_P>(0, *knockout_loc);
                            put<MUTATION_PER_SITE_P>(0, *knockout_loc);
                            put<GERM_MUTATION_PER_SITE_P>(0, *knockout_loc);
                            
                            int cur_update = 0;
                            int update_max = 10000;
                            // and run till the group amasses the right amount of resources
                            while ((get<DIVIDE_REMOTE>(*knockout_loc,0) == 0) &&
                                   (cur_update < update_max)){
                                knockout_loc->update();
                                ++cur_update;
                            }
                            
                            
                            
                            if (knockout_loc->population().size() < 2) {
                                num_uni++;
                                
                                if (cur_update == update_max) {
                                    num_uni_inviable++;
                                    continue;
                                }
                                
                                
                                num_uni_viable++;
                                
                                
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
                                    
                                    typename EA::population_type init_mc;
                                    init_mc.insert(init_mc.end(),knockout_loc2);
                                    
                                    std::swap(metapop.population(), init_mc);
                                    
                                    add_event<mt_gls_propagule>(metapop);
                                    
                                    int max_size = 32;
                                    int max_update = 50000;
                                    cur_update = 0;
                                    
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
                                    
                                    df.write(lod_depth)
                                    .write("uni")
                                    .write(uni_count)
                                    .write(nr)
                                    .write(cur_update)
                                    .write(total_workload)
                                    .write(metapop.size());
                                    df.endl();
                                    
                                    
                                }
                                uni_count++;
                                
                            }
                            
                        }
                        
                    }
                    
                    df2.write(lod_depth)
                    .write(uni_count)
                    .write(num_uni_viable)
                    .write(num_uni_inviable)
                    .write(birth_up)
                    .endl();
                }
                
                
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
    int timepoint = get<ANALYSIS_LOD_TIMEPOINT_TO_ANALYZE>(ea,0);

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
        
    df.write(timepoint)
    .write(i->size())
    .write(num_germ)
    .write(total_workload)
    .write(germ_workload)
    .endl();
    
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
                
                add_ent *= 2;
                if (add_ent >= 1024) {
                    entrench_not_found = false;
                }
                
                
            }// end while
        }
    
    
    LIBEA_ANALYSIS_TOOL(lod_entrench_mc_cost) {
            
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
            int start_rep = get<ANALYSIS_LOD_START_REP>(ea,0);

            int mc_res = get<GROUP_REP_THRESHOLD>(ea);
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
                for (int nr = start_rep; nr < start_rep+num_rep; nr++) {
                    // should define checkpoint + analysis input
                    //ea is the thing loaded from the checkpoint; EA is its type
                    EA metapop; // a new EA
                    
                    metapop.initialize(ea.md());
                    put<GROUP_REP_THRESHOLD>(mc_res, metapop);
                    
                    int new_seed = ea.rng().uniform_integer();
                    put<RNG_SEED>(new_seed, metapop);
                    metapop.reset_rng(new_seed);
                    
                    float start_gen = 0;
                    typename EA::population_type init_mc;
                    for (int j=0; j<meta_size; ++j){
                        typename EA::individual_ptr_type control_mc = metapop.make_individual(*i->traits().founder());
                        control_mc->initialize(metapop.md());
                        put<GROUP_REP_THRESHOLD>(mc_res, *control_mc);
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
                            if (metapop.current_update() > 1000){
                                if (mean_size < 2) {
                                    exit_mean_size++;
                                } else  {
                                    exit_mean_size = 0;
                            }}
                            
                            
                           
                            
                            df.write(mc_res)
                            .write(nr)
                            .write(metapop.current_update())
                            .write(organism_size/metapop.size())
                            .write(num_germ/metapop.size())
                            .write(mean_gen)
                            .write(mean_gen_diff)
                            .write(total_workload/organism_size)
                            .write(germ_workload/num_germ)
                            .endl();
                            
                            if ((exit_mean_size >= 5) ||
                                (mean_gen_diff > 100))  {
                                int reverted = 0;
                                
                                if ((exit_mean_size >= 5) &&
                                    (mean_gen_diff > 10)){
                                    revert_count += 1;
                                    reverted = 1;
                                }
                                exit = true;
                                df2.write(mc_res)
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
                                df2.write(mc_res)
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
                
                mc_res -= 25;
                if (mc_res < 0) {
                    entrench_not_found = false;
                }
                
                
            }// end while
        }
    
    
    
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
        
    }
}


#endif
