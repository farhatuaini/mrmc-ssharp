/**
*	WARNING: Do Not Remove This Section
*
*       $LastChangedRevision: 415 $
*       $LastChangedDate: 2010-12-18 17:21:05 +0100 (Sa, 18. Dez 2010) $
*       $LastChangedBy: davidjansen $
*
*	MRMC is a model checker for discrete-time and continuous-time Markov
*	reward models. It supports reward extensions of PCTL and CSL (PRCTL
*	and CSRL), and allows for the automated verification of properties
*	concerning long-run and instantaneous rewards as well as cumulative
*	rewards.
*
*	Copyright (C) The University of Twente, 2004-2008.
*	Copyright (C) RWTH Aachen, 2008-2009.
*	Authors: Maneesh Khattri, Ivan Zapreev
*
*       Copyright (C) (this file) Saarland University, 2007
*       Authors: Moritz Hahn
*
*	This program is free software; you can redistribute it and/or
*	modify it under the terms of the GNU General Public License
*	as published by the Free Software Foundation; either version 2
*	of the License, or (at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program; if not, write to the Free Software
*	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*
*	Main contact:
*		Lehrstuhl für Informatik 2, RWTH Aachen University
*		Ahornstrasse 55, 52074 Aachen, Germany
*		E-mail: info@mrmc-tool.org
*
*       Old contact:
*		Formal Methods and Tools Group, University of Twente,
*		P.O. Box 217, 7500 AE Enschede, The Netherlands,
*		Phone: +31 53 4893767, Fax: +31 53 4893247,
*		E-mail: mrmc@cs.utwente.nl
*
*       Our contacts (Saarland University).
*               Sven Johr, Dependable Systems and Software
*               Stuhlsatzenhausweg 45
*               66123 Saarbruecken
*               Phone: +49 681 302 5607
*               Fax: +49 681 302 5636
*               Location:   Bldg. E1 3, Room 536
*               E-Mail: johr@cs.uni-sb.de
*
*       Source description: Reads an MDP
*/


#ifndef READ_MDPI_FILE_H
#define READ_MDPI_FILE_H

#include "mdp_sparse.h"

#define MAX_LINE_LENGTH 1024

/**
* Returns the filename of the current MDPI model.
*
* @return filename of current MDPI model
*/
extern const char * mdpi_get_model_name(void);

/**
* Reads CTMDPI from file @a filename.
*
* @param filename file to read CTMDPI from
* @return CTMDPI read from file
*/
extern NDSparseMatrix *read_NDSparseMatrix_file(const char*);

#endif
