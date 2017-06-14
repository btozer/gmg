"""
SAVE CURRENT FaT-FM DISPLAY AS A EDITABLE FIGURE
"""

###########################%
#% SETUP FIGURE PARAMETERS#%
###########################%

#% IMPORT PYTHON MODULES
import numpy as np, copy
import pylab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
import os

def draw_lines(file_name, ax3):

    #% READ IN DATA
    lines = np.genfromtxt(file_name, dtype=str, autostrip=True, delimiter=' ', filling_values='>')

    print "rays loaded"
    print rays
    #% INIT LISTS
    plotx_ray_list = [[]]
    ploty_ray_list = [[]]

    #% EXTRACT LINES
    for i in range(0,len(rays)):
        c=0
        if rays[i, 0] != ">":
            plotx_ray_list[c].append(rays[i,0])
            ploty_ray_list[c].append(rays[i,1])
        elif rays[i, 0] == ">":
            plotx_ray_list.append([])
            ploty_ray_list.append([])
            c=c+1
    #% PLOT LINES
    print "!!!!!!!!!!!!!!"
    print ploty_ray_list[0]
    print "!!!!!!!!!!!!!!"
    for l in range(0, 1):
        ax3.plot(plotx_ray_list[l], ploty_ray_list[l], color='r', linewidth=.08)

########################################################################################################################
########################################################################################################################

def plot_fig(file_path, area, xp, obs_grav, calc_grav, obs_mag, calc_mag, layer_count, layer_lock_list, plotx_list, ploty_list,
             densities, absolute_densitites, reference_densities, segy_plot_list, well_list, well_name_list, t_canvas,
             d_canvas, nt_canvas, model_aspect, use_tight_layout, poly_alpha, fs, ms, lw, font_type, layer_colors,
             draw_polygons, draw_layers, draw_floating_layers, draw_colorbar, draw_xy_data, xy_size, xy_color, colorbar_x, colorbar_y,
             colorbar_size_x, colorbar_size_y, layer_line_width, layer_alpha, grav_rms_value, mag_rms_value, grav_y_min,
             grav_y_max, xy_list_save, draw_wells, wells, well_fs, well_line_width):

    print "starting plot_fig function"
    print "self.segy_plot_list = %s" % segy_plot_list
    print "save_path = %s" % file_path
    ######################################################
    #% SET DEFAULT PLOTTING PARAMS                      %#
    ######################################################
    from matplotlib import rcParams
    #%FIGURE NAME
    fig_name='test'
    #% FIGURE FONT SIZE
    fs = fs
    #% FONT TYPE
    #plt.rc('font', family='serif')
    plt.rc('font', family=font_type, size=fs)
    #% MARKER SIZE
    ms = ms
    #% LINE WIDTH
    lw = lw
    #% AXIS TICK FONT SIZE
    plt.rc('xtick', labelsize=8.0)
    plt.rc('ytick', labelsize=8.0)
    #% AXIS TICK POSITIONS
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    #% SET FONT TYPE AS ONE WHICH IS EDITABLE BY ILLUSTRATOR/INKSCAPE ETC
    rcParams['pdf.fonttype'] = 42
    #% GRID LINE STYLE
    #plt.rc('grid', c='0.5', ls='-', lw=1)
    #from matplotlib.patches import Polygon
    #import custom_hatch_patterns

    '''#%dir containing borehole icon'''
    borehole_dir = os.path.dirname(__file__)+'/icons/'

    #% NUMBER OF ROWS IN PLOT
    if t_canvas is True and d_canvas is True and nt_canvas is True:
        num_rows=12 #% determines axes sizings
    if t_canvas is False and d_canvas is True and nt_canvas is True:
        num_rows=11 #% determines axes sizings
    if t_canvas is True and d_canvas is False and nt_canvas is True:
        num_rows=11 #% determines axes sizings
    if t_canvas is True and d_canvas is True and nt_canvas is False:
        num_rows=11 #% determines axes sizings
    if t_canvas is False and d_canvas is False and nt_canvas is True:
        num_rows=10 #% determines axes sizings
    if t_canvas is True and d_canvas is False and nt_canvas is False:
        num_rows=10 #% determines axes sizings
    if t_canvas is False and d_canvas is True and nt_canvas is False:
        num_rows=10 #% determines axes sizings
    if t_canvas is False and d_canvas is False and nt_canvas is False:
        num_rows=9 #% determines axes sizings

    print "!!!!!!!!!!!!!!!!!!"
    print obs_mag
    print calc_mag
    print "!!!!!!!!!!!!!!!!!!"

    #% DETERMINE AXES
    if obs_grav != []:
        y_start_grav = np.append(obs_grav[:, 1], calc_grav).min()-\
                   abs(np.append(obs_grav[:, 1], calc_grav).min()/20)
        y_end_grav = np.append(obs_grav[:, 1], calc_grav).max()+\
                 (np.append(obs_grav[:, 1], calc_grav).max()/20)


    if obs_mag != []:
        y_start_mag = np.append(obs_mag[:, 1], calc_mag).min()-\
                  abs(np.append(obs_mag[:, 1], calc_mag).min()/20)
        y_end_mag = np.append(obs_mag[:, 1], calc_mag).max()+\
                abs(np.append(obs_mag[:, 1], calc_mag).max()/20)

    x_start_model, x_end_model, y_end_model, y_start_model = np.array(area)
    #print "y_start_grav = %s" % y_start_grav
    #print "y_end_grav = %s" % y_end_grav

    #######################################################
    #% START FIGURE PLOTTING                             #%
    #######################################################
    fig = plt.figure(dpi=600, facecolor='w', edgecolor='k', figsize=(11.69, 8.27 ))
    row_counter = 0

    ####################################################################################################################
    ####################################################################################################################

    '#% 1. PLOT TOPOGRAPHY CANVAS'
    if t_canvas == True:
        ax1 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=1, colspan=1)
        #ax1.set_aspect(model_aspect)
        row_counter=row_counter+1

        '#% AXIS OPTIONS'
        plt.ylabel('Topography. (km)', fontsize=fs)
        plt.xlim(x_start_model, x_end_model)
        plt.ylim(y_start_grav, y_end_grav)
        #all_sub.set_alpha(1.0)
        #all_sub.set_alpha(1.0)
        #ax1.set_yticks([])
        ax1.set_xticks([])
        #% PLOT TOPO DATA
        #ax1.scatter(obs_topo[:, 0], obs_topo[:, 1], marker='o', color='b', s=ms)
        #ax1.plot(xp*0.001, calc_grav, color='r', linewidth=lw)

    '#% 2. PLOT GRAVITY CANVAS'
    if d_canvas == True:
        ax2 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=1, colspan=1)
        #ax2.set_aspect(model_aspect)
        row_counter=row_counter+1
        #plt.grid(zorder=2)
        #uncon_alpha=0.2
        #plt.axvspan(400., 405., facecolor='#B6B6B4', zorder=1, alpha=uncon_alpha) #% Caledonian Orogeney

        '#% AXIS OPTIONS'
        #ax.xaxis.set_label_position('top')
        #plt.xlabel('Distance (km)')
        #plt.grid()
        plt.ylabel('Bouguer Anomaly (mGal)', fontsize=fs, labelpad=-1)
        plt.xlim(x_start_model, x_end_model)
        plt.ylim(grav_y_min, grav_y_max)
        #all_sub.set_alpha(1.0)
        #all_sub.set_alpha(1.0)
        #ax2.set_aspect(obs_aspect_ratio)
        #ax2.set_yticks([])
        #ax2.set_xticks([])
        ax2.set_xticklabels([])
        #fontsize=5

        '#% PLOT GRAVITY DATA'
        if obs_grav != []:
            ax2.scatter(obs_grav[:, 0], obs_grav[:, 1], marker='o', color='b', s=ms)
            #% SET AXIS POSTIONS
            ax2.spines['right'].set_color('none')
            ax2.spines['bottom'].set_position('center')
            ax2.spines['top'].set_position('center')
            ax2.tick_params(axis='x', which='both', labelbottom='off', labeltop='off')
            ax2.tick_params(axis='y', which='both', left='on', right='off', labelright='off')

        if calc_grav is not None:
            ax2.plot(xp*0.001, calc_grav, color='r', linewidth=lw)

        '#% RMS value'
        ax2.annotate('RMS misfit: '+str(grav_rms_value), xy=(10, 60), xytext=(10, 60),
                     fontsize=fs, horizontalalignment='left', clip_on=False)
        if grav_rms_value != 0:
            pass

    '#% 3. PLOT MAGNETIC CANVAS'
    if nt_canvas == True:
        ax3 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=1, colspan=1)
        #ax3.set_aspect(model_aspect)
        row_counter=row_counter+1
        #plt.grid(zorder=2)
        #% PLOT PARNAIBA UNCONFORMITIES
        #uncon_alpha=0.2
        #plt.axvspan(400., 405., facecolor='#B6B6B4', zorder=1, alpha=uncon_alpha) #% Caledonian Orogeney

        '#% AXIS OPTIONS'
        #ax.xaxis.set_label_position('top')
        #plt.xlabel('Distance (km)')
        #plt.grid()
        plt.ylabel('Magnetic \n Anom. (nT)', fontsize=fs)
        plt.xlim(x_start_model, x_end_model)
        plt.ylim(y_start_mag, y_end_mag)
        #all_sub.set_alpha(1.0)
        #all_sub.set_alpha(1.0)
        #ax1.set_yticks([])
        ax3.set_xticks([])
        #fontsize=5

        '#% PLOT MAGNETIC DATA'
        if obs_mag != []:
            ax3.scatter(obs_mag[:, 0], obs_mag[:, 1], marker='o', color='b', s=ms)
        if calc_mag is not None:
            ax3.plot(xp*0.001, calc_mag, color='r', linewidth=lw)

    ###################################################################################################################
    ###################################################################################################################
    print "FORWARD MODEL PANEL"

    '#% 3. FORWARD MODEL PANEL'
    ax4 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=6, colspan=1)
    ax4.set_aspect(model_aspect)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off

    #plt.axvspan(439., 425., facecolor='#798989', alpha=0.5, hatch=mst_hatch, edgecolor='black') #% Tiangua
    #plt.text(439.-((439.-425.)/2.), 45, 'Tia.', fontsize=fontsize_smaller, ha='center', va='center')

    '#% PLOT SEISMIC DATA'
    for s in range(0, len(segy_plot_list)):
        if segy_plot_list[s] != []:
            ax4.add_image(copy.copy(segy_plot_list[s]))

    if draw_wells == True:
        print "DOING WELL"
        '#% PLOT WELL DATA'
        well_textsize = 2.0
        well_lw = 2.0
        for w in range(0, len(well_list)):
            well_data = np.array(well_list[w])
            if well_name_list[w] == "None" or well_name_list[w] == []:
                continue
            else:
                if wells[w][0].get_visible() == False:
                    continue
                else:
                    well_data = np.array(well_list[w])
                    print "!!!!!!!!!!!!!!!!!!!!!!!!"
                    print "well_data = %s" % well_data
                    print "!!!!!!!!!!!!!!!!!!!!!!!!"
                    y1 = well_data[0, 1].astype(float)* -1.
                    y2 = well_data[-1, -1].astype(float)+(y1.astype(float))
                    well_x_location = well_data[1, 1]
                    wellx = (well_x_location, well_x_location)
                    welly = (y1, y2)
                    ax4.plot(wellx, welly, linestyle='-', linewidth=well_line_width, color='black')

                    ax4.annotate(well_name_list[w], xy=(well_x_location, -0.5),
                                                                   xytext=(well_x_location, y1-0.1),
                                                                   fontsize=well_fs, weight='bold',
                                                                   horizontalalignment='center', color='black',
                                                                   bbox=dict(boxstyle="round,pad=.2", fc="0.8",
                                                                             ec='None'), clip_on=True)
                    '''#% PLOT WELL HORIZONS'''
                    for i in range(2, len(well_data)):
                        y = [well_data[i, 1].astype(float)-well_data[0, 1].astype(float),
                             well_data[i, 1].astype(float)-well_data[0, 1].astype(float)]
                        x = [well_data[1, 1].astype(float)-1, well_data[1, 1].astype(float)+1]
                        ax4.plot(x, y, linestyle='-', linewidth=well_line_width, color='black')
                        horizon_y_pos = well_data[i, 1].astype(float)-well_data[0, 1].astype(float) + 0.01
                        horizon = well_data[i, 0].astype(str)

                        '#% ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP'
                        if i % 2 == 0:
                            horizon_x_pos = well_data[1, 1].astype(float) - 3.05
                            ax4.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                    xytext=(horizon_x_pos, horizon_y_pos), fontsize=well_fs,
                                                    weight='bold', horizontalalignment='left', verticalalignment='top',
                                                    color='black', bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'),
                                                                   clip_on=True)
                        else:
                            horizon_x_pos = well_data[1, 1].astype(float) + 3.05
                            ax4.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
                                                    xytext=(horizon_x_pos, horizon_y_pos), fontsize=well_fs,
                                                    weight='bold', horizontalalignment='right', verticalalignment='bottom',
                                                    color='black', bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'),
                                                                   clip_on=True)

    print "PLOT LAYER POLYGON"
    '#% PLOT LAYER POLYGONS'
    if draw_polygons is True:
        for i in range(layer_count, -1, -1):
            #% CREATE POLYGONS
            if layer_lock_list[i] == 0 and i >= 1:
                #%Create polygon
                plotx_polygon = np.append(np.array(plotx_list[i]), np.array(plotx_list[i-1])[::-1])
                ploty_polygon = np.append(np.array(ploty_list[i]), np.array(ploty_list[i-1])[::-1])
            else:
                 plotx_polygon = np.array(plotx_list[i])
                 ploty_polygon = np.array(ploty_list[i])


            '#%DEFINE COLOR MAP'
            colormap = cm.coolwarm
            cNorm  = colors.Normalize(vmin=-0.8, vmax=0.8)
            colormap = cm.ScalarMappable(norm=cNorm, cmap=colormap)

            #% CREATE POLYGON FILL
            if densities[i] != 0 and absolute_densitites == True:
                next_color = colormap.to_rgba(0.001*densities[i] - 0.001*reference_densities[i])
            elif densities[i] != 0:
                next_color = colormap.to_rgba(0.001 * densities[i])
            else:
                next_color = colormap.to_rgba(0.0)

            '#%DRAW'
            ax4.fill(plotx_polygon, ploty_polygon, color=next_color, alpha=poly_alpha, closed=True, ec='none', zorder=1)

    '#% PLOT LAYER LINES'
    if draw_layers is True:
        for i in range(layer_count, -1, -1):
            #% CREATE POLYGONS
            if layer_lock_list[i] == 0 and i >= 1:
                 plotx = np.array(plotx_list[i])
                 ploty = np.array(ploty_list[i])
                 ax4.plot(plotx, ploty, color=layer_colors[i], linewidth=layer_line_width, alpha=layer_alpha, zorder=1)
            else:
                pass

    if draw_floating_layers is True:
        for i in range(layer_count, -1, -1):
            #% CREATE POLYGONS
            if layer_lock_list[i] == 1 and i >= 1:
                 plotx = np.append(plotx_list[i], plotx_list[i][0])
                 ploty = np.append(ploty_list[i], ploty_list[i][0])
                 ax4.plot(plotx, ploty, color=layer_colors[i], linewidth=layer_line_width, alpha=layer_alpha, zorder=1)
            else:
                pass

    '#% DRAW OTHER XY DATA e.g. EARTHQUAKE HYPOCENTERS'
    if draw_xy_data == True:
        print "DRAWING XY DATA"
        for i in xrange(len(xy_list_save)):
            if xy_list_save[i] != []:
                xy = xy_list_save[i]
                ax4.scatter(xy[:, 0], xy[:, 1], marker='o', edgecolors='none', facecolors=xy_color,
                            s=xy_size, gid=i, alpha=1.0, zorder=2)

    '#% PLOT OTHER LINES'
    #ax4.plot([0., 235., 235, 0., 0.], [0., 0., 45., 45., 0.], color='green', linewidth=lw)
    #draw_lines('rays.rayplot', ax3)

    '#% AXIS OPTIONS'
    #ax.xaxis.set_label_position('top')
    plt.xlabel('Distance (km)', fontsize=fs, labelpad=-1)
    plt.ylabel('Depth (Km)', fontsize=fs, labelpad=-1)
    plt.xlim(x_start_model, x_end_model)
    plt.ylim(y_end_model, y_start_model)
    ax4.spines['top'].set_color('none')
    ax4.tick_params(axis='x', which='both', labelbottom='on', labeltop='off')

    ####################################################################################################################
    ####################################################################################################################

    print "SET AXIS DIMENSIONS"
    '#% SET axis dimensions so xaxis is the same as the model plot'
    pos1 = ax4.get_position()
    if t_canvas == True:
        pos2 = ax1.get_position()
        ax1.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])

    if d_canvas == True:
        pos2 = ax2.get_position()
        ax2.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height*2])

    if nt_canvas == True:
        pos2 = ax3.get_position()
        ax3.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])

    ####################################################################################################################
    ####################################################################################################################

    '# PLOT COLOR MAP'
    if draw_colorbar == True:
        colormap      = cm.coolwarm
        cNorm         = colors.Normalize(vmin=-0.8, vmax=0.8)
        C      = cm.ScalarMappable(cmap=colormap, norm=cNorm)
        C._A = []
        c_cax = plt.axes([colorbar_x, colorbar_y, colorbar_size_x, colorbar_size_y])
        cbar = plt.colorbar(C, ticks=[-0.8,-0.4,0,0.4,0.8], orientation="horizontal", cax=c_cax)
        cbar.set_label('Density contrast ($g/cm^{3}$)', fontsize=fs, labelpad=-1)

    ####################################################################################################################
    ####################################################################################################################

    '#% WRITE OUT FIG'
    if use_tight_layout == True:
        # plt.savefig(outdir+fig_name+'.svg', bbox_inches='tight', dpi=720, format='svg')
        # print "SAVED SVG"
        plt.savefig(file_path, bbox_inches='tight', dpi=720, format='pdf')
        print "SAVED PDF"
        # plt.savefig(outdir+fig_name+'.ps', bbox_inches='tight', dpi=720, format='ps')
        # print "SAVED PS"
        # plt.savefig(outdir+fig_name+'.eps', bbox_inches='tight', dpi=720, format='eps')
        # print "SAVED EPS"
    else:
        # plt.savefig(outdir+fig_name+'.svg', bbox_inches='tight', dpi=720, format='svg')
        # print "SAVED SVG"
        plt.savefig(file_path, dpi=720, format='pdf')
        print "SAVED PDF"
        # plt.savefig(outdir+fig_name+'.ps', bbox_inches='tight', dpi=720, format='ps')
        # print "SAVED PS"
        # plt.savefig(outdir+fig_name+'.eps', bbox_inches='tight', dpi=720, format='eps')
        # print "SAVED EPS"

    print "figure saved"
    return
########################################################################################################################
########################################################################################################################
