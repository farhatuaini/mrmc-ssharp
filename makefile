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
#	Source description: The main Markov Reward Model Checker
#				m a k e f i l e.
#

export MRMC_HOME_DIR = $(PWD)

.PHONY: error all lib lint depend clean

error:
	@echo "+-----------------------------------------------------------------+"
	@echo "|                                                                 |"
	@echo "|          Install Markov Reward Model Checker (MRMC)             |"
	@echo "|                                                                 |"
	@echo "|       Usage: make all    - installs MRMC                        |"
	@echo "|              make clean  - cleans *.o and other executables     |"
	@echo "|                                                                 |"
	@echo "|  Make sure the system-specific makefile.def has been edited     |"
	@echo "|  to reflect your system configuration.                          |"
	@echo "|                                                                 |"
	@echo "+-----------------------------------------------------------------+"


bin:
	mkdir $@

lib:
	-mkdir $@
	$(MAKE) -C obj $@

all: bin lib
	$(MAKE) -C obj $@

lint:
	$(MAKE) -C obj $@

depend:
	$(MAKE) -C obj $@

clean:
	$(MAKE) -C obj $@
	rm -rf lib
