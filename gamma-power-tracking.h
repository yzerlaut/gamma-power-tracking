/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This is a template header file for a user modules derived from
 * DefaultGUIModel with a custom GUI.
 */

#include <default_gui_model.h>
#include <math.h>

const int Length_LFP_vector = 2000; 
const int Length_pLFP_vector = 10; 
const double w0 = 6.; // parameter for wavelet extent

class GammaPowerTracking : public DefaultGUIModel
{

  Q_OBJECT


public:

  GammaPowerTracking(void);
  virtual ~GammaPowerTracking(void);

  void execute(void);
  void createGUI(DefaultGUIModel::variable_t*, int);
  void increase_Real_and_Imaginary_of_Morlet_TF(double*, double*, double,
						double, double, double, double,
						int, int);
  void normalization_factor(double*,
			    double, double, double);


protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  void initParameters();
  void initWavelets();

  /* LFP props */
  double LFP_history_vector[Length_LFP_vector];
  double pLFP_history_vector[Length_pLFP_vector];
  double pLFP_history_norm_vector[Length_pLFP_vector];
  double mean_LFP_last;
  double mean_LFP_new;
  double cum_pLFP;

  /* WAVELET */
  double wavelet1_freq;
  double wavelet1_norm;
  /* double wavelet2_freq; */
  int Length_wavelet1;
  int Length_wavelet2;

  /* SMOOTHING */
  double Tsmoothing;
  int ismoothing;

  /* SYSTEM */
  double period;
  long long count;
  double systime;


};
