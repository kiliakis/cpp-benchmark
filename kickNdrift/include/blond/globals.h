/*
 * globals.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_GLOBALS_H_
#define INCLUDES_GLOBALS_H_

#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Beams.h>
#include <blond/beams/Slices.h>

struct API Context {
    static Slices *Slice;
    static GeneralParameters *GP;
    static Beams *Beam;
    static RfParameters *RfP;

    // TODO num of threads is global
    // should it?
    static int n_threads;
};

#endif /* INCLUDES_GLOBALS_H_ */
