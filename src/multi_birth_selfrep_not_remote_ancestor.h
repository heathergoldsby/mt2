/* digital_evolution/ancestors/multi_birth_selfrep_not_remote_ancestor.h
 *
 * This file is part of EALib.
 *
 * Copyright 2012 David B. Knoester, Heather J. Goldsby.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _EA_DIGITAL_EVOLUTION_ANCESTORS_MULTI_BIRTH_SELFREP_NOT_REMOTE_ANCESTOR_H_
#define _EA_DIGITAL_EVOLUTION_ANCESTORS_MULTI_BIRTH_SELFREP_NOT_REMOTE_ANCESTOR_H_

#include <ea/digital_evolution.h>
#include <ea/metadata.h>


namespace ealib {
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct multibirth_selfrep_not_remote_ancestor {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["fixed_input"]; // input
            repr[25] = ea.isa()["fixed_input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
//h_alt_divide            repr[96] = ea.isa()["h_alt_divide"];
            repr[96] =  ea.isa()["h_divide_remote"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };

    /*! Generates a self-replicating ancestor that performs not.
     */
    struct multibirth_selfrep_not_remote_ancestor_unfixed {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["input"]; // input
            repr[25] = ea.isa()["input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            //h_alt_divide            repr[96] = ea.isa()["h_alt_divide"];
            repr[96] =  ea.isa()["h_divide_remote"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };

    
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct multibirth_selfrep_prop1_remote_ancestor {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            /*
            repr[10] = ea.isa()["if_not_birth_0"];
            repr[11] = ea.isa()["become_soma"];
            repr[13] = ea.isa()["if_not_birth_0"];
            repr[14] = ea.isa()["flag_1"];
             */
            repr[10] = ea.isa()["if_not_member_start_propagule"];
            repr[11] = ea.isa()["become_soma"];
            repr[13] = ea.isa()["if_not_member_start_propagule"];
            repr[14] = ea.isa()["flag_1"];
            
            // not
            repr[24] = ea.isa()["fixed_input"]; // input
            repr[25] = ea.isa()["fixed_input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            //h_alt_divide            repr[96] = ea.isa()["h_alt_divide"];
            repr[96] =  ea.isa()["h_divide_remote"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };

    
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct multibirth_selfrep_not_remote_ancestor_mod_resources {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["fixed_input"]; // input
            repr[25] = ea.isa()["fixed_input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output_res_moderated_by_group_size"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            //h_alt_divide            repr[96] = ea.isa()["h_alt_divide"];
            repr[96] =  ea.isa()["h_divide_remote"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };
    
    struct multibirth_selfrep_not_remote_comm_ancestor {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["fixed_input"]; // input
            repr[25] = ea.isa()["fixed_input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["tx_msg_check_task"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            //h_alt_divide            repr[96] = ea.isa()["h_alt_divide"];
            repr[96] =  ea.isa()["h_divide_remote"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };
    
    struct multibirth_selfrep_not_remote_unfixed_ancestor {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["input"]; // input
            repr[25] = ea.isa()["input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            //h_alt_divide            repr[96] = ea.isa()["h_alt_divide"];
            repr[96] =  ea.isa()["h_divide_remote"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };
    
}

#endif
