/*
 * atomiRadii.h
 *
 *  Created on: Jul 5, 2018
 *      Author: cmyers
 *	hhhh
 */

#ifndef ATOMICRADII_H_
#define ATOMICRADII_H_

struct atomicRadii{
	
	static double bohrRadius_ang;

	static double BW_small(int atom)
	{
		switch(atom)
		{
			case 1:  return 1.45;
			case 6:  return 1.5;
			case 7:  return 2.0;
			case 8:  return 2.0;
			//the following are our own radii;
			//these are not in the paper.
			case 15: return 2.5;
			case 16: return 2.5;
			default: return 2.0;
		}
	}

	static double BW_large(int atom)
	{
		switch(atom)
		{
			case 1:  return 1.5;
			case 6:  return 2.0;
			case 7:  return 2.0;
			case 8:  return 2.0;
			//the following are our own radii;
			//these are not in the paper.
			case 15: return 2.5;
			case 16: return 2.5;
			default: return 2.0;
		}
	}

    static double Bondi_radii(int atom)
    {
        switch (atom)
        {
        case 1:  return 1.06;
        case 6:  return 1.53;
        case 7:  return 1.46;
        case 8:  return 1.42;
            //the following are our own radii;
            //these are not in the paper.
        case 15: return 1.86;
        case 16: return 1.80;
        default: return 2.0;
        }
    }
};




#endif /* ATOMICRADII_H_ */
