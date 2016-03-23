#
#	WARNING: Do Not Remove This Section
#
#       $LastChangedRevision: 415 $
#       $LastChangedDate: 2010-12-18 17:21:05 +0100 (Sa, 18. Dez 2010) $
#       $LastChangedBy: davidjansen $
#
#	MRMC is a model checker for discrete-time and continuous-time Markov
#	reward models. It supports reward extensions of PCTL and CSL (PRCTL
#	and CSRL), and allows for the automated verification of properties
#	concerning long-run and instantaneous rewards as well as cumulative
#	rewards.
#
#	Copyright (C) The University of Twente, 2004-2008.
#	Copyright (C) RWTH Aachen, 2008-2009.
#       Authors: Maneesh Khattri, Ivan Zapreev, David N. Jansen
#
#	This program is free software; you can redistribute it and/or
#	modify it under the terms of the GNU General Public License
#	as published by the Free Software Foundation; either version 2
#	of the License, or (at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program; if not, write to the Free Software
#	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#	MA  02110-1301, USA.
#
#	Main contact:
#		Lehrstuhl für Informatik 2, RWTH Aachen University
#		Ahornstrasse 55, 52074 Aachen, Germany
#		E-mail: info@mrmc-tool.org
#
#       Old contact:
#		Formal Methods and Tools Group, University of Twente,
#		P.O. Box 217, 7500 AE Enschede, The Netherlands,
#		Phone: +31 53 4893767, Fax: +31 53 4893247,
#		E-mail: mrmc@cs.utwente.nl
#
#	Source description: The MRMC makefile that does the hard work.
#


# execute "make inonego"


MRMC_HOME_DIR ?= $(CURDIR)

include $(MRMC_HOME_DIR)/makefile.def

SRC_DIR = $(MRMC_HOME_DIR)/src

LIB_SRC =	$(SRC_DIR)/algorithms/bscc.c \
	$(SRC_DIR)/algorithms/foxglynn.c \
	$(SRC_DIR)/algorithms/iterative_solvers.c
LIB_SRC +=	$(SRC_DIR)/algorithms/random_numbers/rand_num_generator.c \
	$(SRC_DIR)/algorithms/random_numbers/rng_app_crypt.c \
	$(SRC_DIR)/algorithms/random_numbers/rng_ciardo.c \
	$(SRC_DIR)/algorithms/random_numbers/rng_gsl.c \
	$(SRC_DIR)/algorithms/random_numbers/rng_prism.c \
	$(SRC_DIR)/algorithms/random_numbers/rng_ymer.c
LIB_SRC +=	$(SRC_DIR)/io/read_impulse_rewards.c \
	$(SRC_DIR)/io/read_lab_file.c \
	$(SRC_DIR)/io/read_rewards.c \
	$(SRC_DIR)/io/read_tra_file.c \
	$(SRC_DIR)/io/read_mdpi_file.c \
	$(SRC_DIR)/io/execute_cmd_script.c \
	$(SRC_DIR)/io/write_res_file.c \
	$(SRC_DIR)/io/token.c
LIB_SRC +=	$(SRC_DIR)/io/parser/core_to_core.c \
	$(SRC_DIR)/io/parser/parser_to_core.c \
	$(SRC_DIR)/io/parser/parser_to_tree.c
LIB_SRC +=	$(SRC_DIR)/lumping/lump.c \
        $(SRC_DIR)/lumping/sort.c \
	$(SRC_DIR)/lumping/partition.c
LIB_SRC +=	$(SRC_DIR)/modelchecking/prctl.c \
	$(SRC_DIR)/modelchecking/simulation_common.c \
	$(SRC_DIR)/modelchecking/simulation_ctmc.c \
	$(SRC_DIR)/modelchecking/simulation_utils.c \
	$(SRC_DIR)/modelchecking/simulation.c \
	$(SRC_DIR)/modelchecking/steady.c \
	$(SRC_DIR)/modelchecking/transient_common.c \
	$(SRC_DIR)/modelchecking/transient_ctmc.c \
	$(SRC_DIR)/modelchecking/transient_ctmrm.c \
	$(SRC_DIR)/modelchecking/transient_dtmc.c \
	$(SRC_DIR)/modelchecking/transient_dtmrm.c \
	$(SRC_DIR)/modelchecking/transient.c \
	$(SRC_DIR)/modelchecking/transient_ctmdpi_hd_uni.c \
	$(SRC_DIR)/modelchecking/transient_ctmdpi_hd_non_uni.c
LIB_SRC +=	$(SRC_DIR)/storage/bitset.c \
	$(SRC_DIR)/storage/kjstorage.c \
	$(SRC_DIR)/storage/label.c \
	$(SRC_DIR)/storage/mdp_labelset.c \
	$(SRC_DIR)/storage/path_graph.c \
	$(SRC_DIR)/storage/sample_vec.c \
	$(SRC_DIR)/storage/sparse.c \
	$(SRC_DIR)/storage/mdp_sparse.c \
	$(SRC_DIR)/storage/stack.c
LIB_SRC +=	$(SRC_DIR)/runtime.c

NONLIB_SRC =	$(SRC_DIR)/mcc.c

SRC	= $(LIB_SRC) $(NONLIB_SRC)

LEX_SRC = $(SRC_DIR)/io/parser/mrmc_tokenizer.l
LEX_TGT = $(MRMC_HOME_DIR)/obj/$(notdir $(LEX_SRC:.l=.c))

YACC_SRC = $(SRC_DIR)/io/parser/mrmc_grammar.y
YACC_TGT = $(MRMC_HOME_DIR)/obj/$(notdir $(YACC_SRC:%.y=%.tab.c))
YACCH = $(YACC_TGT:.c=.h)

LIB_SRC +=	$(LEX_TGT) $(YACC_TGT)

LIB_OBJ = $(addprefix $(MRMC_HOME_DIR)/obj/,$(notdir $(LIB_SRC:.c=.o)))
NONLIB_OBJ = $(addprefix $(MRMC_HOME_DIR)/obj/,$(notdir $(NONLIB_SRC:.c=.o)))
OBJ = $(LIB_OBJ) $(NONLIB_OBJ)
DEP = $(OBJ:.o=.d)

.PHONY: all lib depend lint clean

all: $(EXEC) ;

$(EXEC): $(NONLIB_OBJ) $(LIB_A)
	@echo === ld ===
	@cd $(MRMC_HOME_DIR); \
	$(CC) $(CFLAGS) -o $@ $(^:$(MRMC_HOME_DIR)/%=%) $(LDFLAGS)

lib: $(LIB_A) ;

$(LIB_A): $(LIB_OBJ)
	@echo === $(AR) $(@:$(MRMC_HOME_DIR)/%=%) ===
	@cd $(MRMC_HOME_DIR); \
	$(AR) $(ARFLAGS) crv $@ $(?:$(MRMC_HOME_DIR)/%=%); \
	$(RANLIB) $(RANLIBFLAGS) $@

depend: $(DEP) ;

clean:
	rm -f $(DEP) $(OBJ) $(LIB_A) $(EXEC) $(LEX_TGT) $(YACC_TGT) $(YACCH)

$(YACC_TGT): $(YACC_SRC)
	@echo === $(YACC) $(<:$(MRMC_HOME_DIR)/%=%) ===
	@echo $(YACC) -d $(YFLAGS) -o $(@:$(MRMC_HOME_DIR)/%=%) $(<:$(MRMC_HOME_DIR)/%=%)
	@cd $(MRMC_HOME_DIR); \
	$(YACC) -d $(YFLAGS) -o $(@:$(MRMC_HOME_DIR)/%=%) $(<:$(MRMC_HOME_DIR)/%=%)

$(YACCH): $(YACC_TGT) ;

# We cannot rely on implicit rules because the directories are different

$(LEX_TGT): $(LEX_SRC)
	@echo === $(LEX) $(<:$(MRMC_HOME_DIR)/%=%) ===
	@echo $(LEX) -t $(LFLAGS) $(<:$(MRMC_HOME_DIR)/%=%) \> $@
	@cd $(MRMC_HOME_DIR); \
	rm -f $@; \
	$(LEX) -t $(LFLAGS) $(<:$(MRMC_HOME_DIR)/%=%) > $@
	
define createrules
$(MRMC_HOME_DIR)/obj/$(notdir $(1:.c=.d)) : $(1)
	@echo === makedepend $$(<:$(MRMC_HOME_DIR)/%=%) ===
	@set -e; rm -f $$@; \
	$(CC) -MM $(CPPFLAGS) $(CFLAGS) -MF $$@.$$$$ $$<; \
	sed 's,$(notdir $(1:.c=\.o))[ :]*,$$(@:.d=.o) $$@ : ,' < $$@.$$$$ > $$@; \
	rm -f $$@.$$$$

$(MRMC_HOME_DIR)/obj/$(notdir $(1:.c=.o)) : $(1)
	@echo === $(CC) $$(<:$(MRMC_HOME_DIR)/%=%) ===
	@cd $(MRMC_HOME_DIR); \
	$(CC) -c $(CPPFLAGS:-I$(MRMC_HOME_DIR)/%=-I%) $(CFLAGS) -o $$@ $$(<:$(MRMC_HOME_DIR)/%=%)

lint lint-$(notdir $(1:.c=)) :: $(1)
	@echo === $(LINT) $$(<:$(MRMC_HOME_DIR)/%=%) ===
	@cd $(MRMC_HOME_DIR); \
	$(LINT) $(LINTFLAGS:-I$(MRMC_HOME_DIR)/%=-I%) $$(<:$(MRMC_HOME_DIR)/%=%)
endef


inonego:
	bison -d -o obj/mrmc_grammar.tab.c src/io/parser/mrmc_grammar.y
	flex -t src/io/parser/mrmc_tokenizer.l > /home/Johannes/mrmc/obj/../obj/mrmc_tokenizer.c
	$(CC) $(CPPFLAGS:-I$(MRMC_HOME_DIR)/%=-I%) $(CFLAGS) $(SRC) -lgsl -lgslcblas -lm -L$(GSL_HOME)/lib -o mrmc

