"""
SAVE GMG MODEL AS A EDITABLE VECTOR GRAPHICS FIGURE.

CALLED FROM gmg draw_model() FUNC
"""

import numpy as np, copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import os
from matplotlib import rcParams

# def draw_line(file_name, ax3):
#     """DRAW EXTERNAL XY LINES ON FIGURE"""
#     # READ IN DATA
#     line = np.genfromtxt(file_name, dtype=str, autostrip=True, delimiter=' ', filling_values='>')
#
#     # INIT LISTS
#     plotx_ray_list = [[]]
#     ploty_ray_list = [[]]
#
#     # EXTRACT LINES. ">" IS USED TO INDICATE A NEW LINE SEGMENT
#     r = 0
#     for i in range(0, len(line)):
#         if line[i, 0] != ">":
#             plotx_ray_list[r].append(line[i, 0])
#             ploty_ray_list[r].append(line[i, 1])
#         elif line[i, 0] == ">":
#             # START NEW LINE
#             plotx_ray_list.append([])
#             ploty_ray_list.append([])
#             r += 1
#
#     # PLOT LINES
#     for l in range(0, 1):
#         ax3.plot(plotx_ray_list[l], ploty_ray_list[l], color='r', linewidth=.08)


def plot_fig(file_path, file_type, use_tight_layout, fs, aspect_ratio, ps, calc_line_width, font_type,
             t_frame_min, t_frame_max, grav_frame_min, grav_frame_max, mag_frame_min, mag_frame_max,
             draw_polygons, polygon_alpha, draw_fixed_layers, layer_line_width, draw_floating_layers, layer_line_alpha,
             draw_colorbar, colorbar_x, colorbar_y, colorbar_size_x, colorbar_size_y, draw_xy_data, xy_size, xy_color,
             draw_wells, well_fs, well_line_width, draw_faults, faults_lw, layer_list, fault_list,
             observed_topography_list, observed_gravity_list, observed_magnetic_list, outcrop_data_list, well_data_list,
             segy_data_list, t_frame, g_frame, m_frame, predicted_gravity, gravity_rms_value,
             predicted_nt, magnetic_rms_value, area, xp):
    """
    CREATE A FIGURE FROM THE RRENT GMG MODEL USING USER SPECIFIED PLOTTING OPTIONS AND SAVE TO AN EXTERNAL FILE
    """

    # ------------------------------------------------------------------------------------------------------------------
    # SET DEFAULT PLOTTING PARAMS

    # FIGURE FONT SIZE
    fs = fs

    print(str(font_type))

    # FONT TYPE
    plt.rc('font', family=str(font_type), size=fs)

    # MARKER SIZE
    ps = ps

    # LINE WIDTH
    lw = calc_line_width

    # AXIS TICK FONT SIZE
    plt.rc('xtick', labelsize=8.0)
    plt.rc('ytick', labelsize=8.0)

    # AXIS TICK POSITIONS
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'

    # SET FONT TYPE AS ONE WHICH IS EDITABLE BY VECTOR PROGRAMS ILLUSTRATOR/INKSCAPE ETC
    rcParams['pdf.fonttype'] = 42

    # GRID LINE STYLE
    plt.rc('grid', c='0.5', ls='-', lw=1)


    # DIR CONTAINING BOREHOLE ICON
    borehole_dir = os.path.dirname(os.path.abspath(__file__))+'../docs/icons/'

    # NUMBER OF ROWS IN PLOT
    if t_frame.get_visible() is True and g_frame.get_visible() is True and m_frame.get_visible() is True:
        num_rows = 12  # DETERMINES AXES SIZING
    if t_frame.get_visible() is False and g_frame.get_visible() is True and m_frame.get_visible() is True:
        num_rows = 11  # DETERMINES AXES SIZING
    if t_frame.get_visible() is True and g_frame.get_visible() is False and m_frame.get_visible() is True:
        num_rows = 11  # DETERMINES AXES SIZING
    if t_frame.get_visible() is True and g_frame.get_visible() is True and m_frame.get_visible() is False:
        num_rows = 11  # DETERMINES AXES SIZING
    if t_frame.get_visible() is False and g_frame.get_visible() is False and m_frame.get_visible() is True:
        num_rows = 10  # DETERMINES AXES SIZING
    if t_frame.get_visible() is True and g_frame.get_visible() is False and m_frame.get_visible() is False:
        num_rows = 10  # DETERMINES AXES SIZING
    if t_frame.get_visible() is False and g_frame.get_visible() is True and m_frame.get_visible() is False:
        num_rows = 10  # DETERMINES AXES SIZING
    if t_frame.get_visible() is False and g_frame.get_visible() is False and m_frame.get_visible() is False:
        num_rows = 9  # DETERMINES AXES SIZING

    # ------------------------------------------------------------------------------------------------------------------

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # START FIGURE PLOTTING

    fig = plt.figure(dpi=720, facecolor='w', edgecolor='k', figsize=(11.69, 8.27))
    row_counter = 0

    # SET THE MODEL AXES
    x_start_model, x_end_model, y_end_model, y_start_model = np.array(area)/1000

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT TOPOGRAPHY PANEL
    if t_frame.get_visible() is True:
        ax1 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=1, colspan=1)
        row_counter += 1

        # AXIS OPTIONS
        plt.ylabel('(km)', fontsize=fs)
        plt.xlim(x_start_model, x_end_model)
        plt.ylim(t_frame_min, t_frame_max)
        ax1.set_xticks([])

        # PLOT TOPO DATA
        if observed_topography_list is not None:
            for i in range(0, len(observed_topography_list)):
                if observed_topography_list[i] is not None:
                    ax1.scatter(observed_topography_list[i].data[:, 0], observed_topography_list[i].data[:, 1],
                                marker='o', color=observed_topography_list[i].color, s=ps)

                # SET AXIS POSTIONS
                ax1.spines['right'].set_color('none')
                ax1.spines['bottom'].set_position('center')
                ax1.spines['top'].set_position('center')
                ax1.tick_params(axis='x', which='both', labelbottom='off', labeltop='off')
                ax1.tick_params(axis='y', which='both', left='on', right='off', labelright='off')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT GRAVITY PANEL
    if g_frame.get_visible() is True:
        ax2 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=1, colspan=1)
        row_counter += 1

        # AXIS OPTIONS
        plt.ylabel('mGal', fontsize=fs, labelpad=-1)
        plt.xlim(x_start_model, x_end_model)
        plt.ylim(grav_frame_min, grav_frame_max)
        ax2.set_xticklabels([])

        # PLOT OBSERVED GRAVITY DATA
        if observed_gravity_list is not None:
            for i in range(0, len(observed_gravity_list)):
                if observed_gravity_list[i] is not None:
                    ax2.scatter(observed_gravity_list[i].data[:, 0], observed_gravity_list[i].data[:, 1], marker='o',
                        color=observed_gravity_list[i].color, s=ps)

        # PLOT PREDICTED GRAVITY ANOMALY
        if predicted_gravity is not None:
            ax2.plot(xp*0.001, predicted_gravity, color='r', alpha=0.25, linewidth=calc_line_width)

        # PLOT RMS VALUE
        if gravity_rms_value is not None:
            ax2.annotate('RMS misfit: ' + str(gravity_rms_value), xy=(10, 60), xytext=(10, 60),
                         fontsize=fs, horizontalalignment='left', clip_on=False)

        # SET AXIS POSITIONS
        ax2.spines['right'].set_color('none')
        ax2.spines['bottom'].set_position('center')
        ax2.spines['top'].set_position('center')
        ax2.tick_params(axis='x', which='both', labelbottom='off', labeltop='off')
        ax2.tick_params(axis='y', which='both', left='on', right='off', labelright='off')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT MAGNETIC PANEL
    if m_frame.get_visible() is True:
        ax3 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=1, colspan=1)
        row_counter += 1

        # AXIS OPTIONS
        plt.ylabel('nT', fontsize=fs)
        plt.xlim(x_start_model, x_end_model)
        plt.ylim(mag_frame_min, mag_frame_max)
        ax3.set_xticks([])

        # PLOT OBSERVED MAGNETIC DATA
        if observed_magnetic_list is not None:
            for i in range(0, len(observed_magnetic_list)):
                if observed_magnetic_list[i] is not None:
                    ax2.scatter(observed_magnetic_list[i].data[:, 0], observed_magnetic_list[i].data[:, 0], marker='o',
                        color=observed_magnetic_list[i].color, s=ps)

        # PLOT PREDICTED MAGNETIC ANOMALY
        if predicted_nt is not None:
            ax2.plot(xp * 0.001, predicted_nt, color='g', linewidth=calc_line_width)

        # PLOT RMS VALUE
        if magnetic_rms_value is not None:
            ax2.annotate('RMS misfit: ' + str(magnetic_rms_value), xy=(10, 60), xytext=(10, 60),
                         fontsize=fs, horizontalalignment='left', clip_on=False)

        # SET AXIS POSITIONS
        ax3.spines['right'].set_color('none')
        ax3.spines['bottom'].set_position('center')
        ax3.spines['top'].set_position('center')
        ax3.tick_params(axis='x', which='both', labelbottom='off', labeltop='off')
        ax3.tick_params(axis='y', which='both', left='on', right='off', labelright='off')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT MODEL PANEL
    ax4 = plt.subplot2grid((num_rows, 1), (row_counter, 0), rowspan=6, colspan=1)
    ax4.set_aspect(aspect_ratio)

    plt.tick_params(axis='x',  # CHANGES APPLY TO THE X-AXIS
                    which='both',  # BOTH MAJOR AND MINOR TICKS ARE AFFECTED
                    bottom='on',  # TICKS ALONG THE BOTTOM EDGE ARE OFF
                    top='off',   # TICKS ALONG THE TOP EDGE ARE OFF
                    labelbottom='on')  # LABELS ALONG THE BOTTOM EDGE ARE OFF

    # AXIS OPTIONS
    plt.xlabel('Distance (km)', fontsize=fs, labelpad=-1)
    plt.ylabel('Depth (Km)', fontsize=fs, labelpad=-1)
    plt.xlim(x_start_model, x_end_model)
    plt.ylim(y_end_model, y_start_model)
    plt.gca().invert_yaxis()
    ax4.spines['top'].set_color('none')
    ax4.tick_params(axis='x', which='both', labelbottom='on', labeltop='off')

    # ******************************************************************************************************************
    # # PLOT SEISMIC DATA
    for s in range(0, len(segy_data_list)):
        if segy_data_list[s].pml_actor.get_visible() is True:
            ax4.add_image(copy.copy(segy_data_list[s]))
    # ******************************************************************************************************************
    # draw_polygons, polygon_alpha, draw_fixed_layers, layer_line_width, draw_floating_layers,
    # ******************************************************************************************************************
    # # PLOT LAYER POLYGONS
    print(draw_polygons)
    print(len(layer_list))
    if draw_polygons is True:
        print("draw_polygons is True")
        for i in range(0, len(layer_list)):
            print(("drawing polygon %s") % i)
            # GET ACTOR ATTRIBUTES
            x = layer_list[i].polygon_mpl_actor[0].get_xy()[:,0]
            y = layer_list[i].polygon_mpl_actor[0].get_xy()[:,1]
            fc =  layer_list[i].polygon_mpl_actor[0].get_fc()
            fa =  layer_list[i].polygon_mpl_actor[0].get_alpha()
            lc = layer_list[i].color

            print(x)
            print(y)
            print(fc)
            print(fa)
            print(lc)

            # DRAW ACTOR LINE
            ax4.plot(x, y, color=lc, linewidth=layer_line_width, alpha=layer_line_alpha, zorder=1)

            # DRAW ACTOR FILL
            ax4.fill(x, y, color=fc, alpha=fa, closed=True, ec='none', zorder=1)
    # ******************************************************************************************************************

    # ******************************************************************************************************************
    # if draw_wells:
    #     # PLOT WELL DATA
    #     for w in range(0, len(well_list)):
    #         if well_name_list[w] == "None" or well_name_list[w] == []:
    #             continue
    #         else:
    #             if not wells[w][0].get_visible():
    #                 continue
    #             else:
    #                 well_data = np.array(well_list[w])
    #                 y1 = well_data[0, 1].astype(float)* -1.
    #                 y2 = well_data[-1, -1].astype(float)+(y1.astype(float))
    #                 well_x_location = well_data[1, 1]
    #                 wellx = (well_x_location, well_x_location)
    #                 welly = (y1, y2)
    #                 ax4.plot(wellx, welly, linestyle='-', linewidth=well_line_width, color='black')
    #
    #                 ax4.annotate(well_name_list[w], xy=(well_x_location, -0.5),
    #                                                                xytext=(well_x_location, y1-0.1),
    #                                                                fontsize=well_fs, weight='bold',
    #                                                                horizontalalignment='center', color='black',
    #                                                                bbox=dict(boxstyle="round,pad=.2", fc="0.8",
    #                                                                          ec='None'), clip_on=True)
    #                 # PLOT WELL HORIZONS
    #                 for i in range(2, len(well_data)):
    #                     y = [well_data[i, 1].astype(float)-well_data[0, 1].astype(float),
    #                          well_data[i, 1].astype(float)-well_data[0, 1].astype(float)]
    #                     x = [well_data[1, 1].astype(float)-1, well_data[1, 1].astype(float)+1]
    #                     ax4.plot(x, y, linestyle='-', linewidth=well_line_width, color='black')
    #                     horizon_y_pos = well_data[i, 1].astype(float)-well_data[0, 1].astype(float) + 0.01
    #                     horizon = well_data[i, 0].astype(str)
    #
    #                     # ALTERNATE POSITION OF ODDs/EVENs TO TRY AND AVOID OVERLAP
    #                     if i % 2 == 0:
    #                         horizon_x_pos = well_data[1, 1].astype(float) - 3.05
    #                         ax4.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
    #                                               xytext=(horizon_x_pos, horizon_y_pos), fontsize=well_fs,
    #                                               weight='bold', horizontalalignment='left', verticalalignment='top',
    #                                               color='black', bbox=dict(boxstyle="round,pad=.4", fc="0.8",
    #                                               ec='None'), clip_on=True)
    #                     else:
    #                         horizon_x_pos = well_data[1, 1].astype(float) + 3.05
    #                         ax4.annotate(horizon, xy=(horizon_x_pos, horizon_y_pos),
    #                                               xytext=(horizon_x_pos, horizon_y_pos), fontsize=well_fs,
    #                                               weight='bold', horizontalalignment='right',
    #                                               verticalalignment='bottom', color='black',
    #                                               bbox=dict(boxstyle="round,pad=.4", fc="0.8", ec='None'), clip_on=True)
    # ******************************************************************************************************************



    # ******************************************************************************************************************
    # if draw_faults is True:
    #     # PLOT FAULTS
    #     for i in range(0, len(faults)-1):
    #         fault = faults[i][0]
    #         # CHECK IF FAULT IS SET AS VISIBLE; IF YES, THEN PLOT
    #         if fault.get_visible() == True:
    #             x = fault.get_xdata()
    #             y = fault.get_ydata()
    #             ax4.plot(x, y, color='k', linewidth=0.5, zorder=1, alpha=1.0)
    #         else:
    #             continue
    # ******************************************************************************************************************

    # ******************************************************************************************************************
    # # DRAW OTHER XY DATA e.g. EARTHQUAKE HYPOCENTERS
    # if draw_xy_data:
    #     for i in range(0, len(xy_list)):
    #         if xy_list[i]:
    #             xy = xy_list[i]
    #             ax4.scatter(xy[:, 0], xy[:, 1], marker='o', edgecolors='none', facecolors=xy_color,
    #                         s=xy_size, gid=i, alpha=1.0, zorder=2)
    # ******************************************************************************************************************

    # ******************************************************************************************************************
    # # PLOT EXTERNAL LINES
    # # for i in range(0, len(external_lines)):
    # #     # READ FILE PATH OF LINE'
    # #     line_file_path = external_lines[i]
    # #     #  DRAW LINE ON FIGURE'
    # #     draw_line(line_file_path, ax3)
    # ******************************************************************************************************************


    # ******************************************************************************************************************
    # # SET AXIS DIMENSIONS SO X AXIS IS THE SAME AS THE MODEL PLOT
    # pos1 = ax4.get_position()
    # if self.t_frame and obs_topo:
    #     pos2 = ax1.get_position()
    #     ax1.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    # if self.g_frame and obs_grav:
    #     pos2 = ax2.get_position()
    #     ax2.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height*2])
    # ifself.m_frame and obs_mag:
    #     pos2 = ax3.get_position()
    #     ax3.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    # ******************************************************************************************************************

    # ******************************************************************************************************************
    # # PLOT COLOR MAP
    # if draw_colorbar:
    #     colormap = cm.coolwarm
    #     cnorm = colors.Normalize(vmin=-0.8, vmax=0.8)
    #     C = cm.ScalarMappable(cmap=colormap, norm=cnorm)
    #     C._A = []
    #     c_cax = plt.axes([colorbar_x, colorbar_y, colorbar_size_x, colorbar_size_y])
    #     cbar = plt.colorbar(C, ticks=[-0.8, -0.4, 0.0, 0.4, 0.8], orientation="horizontal", cax=c_cax)
    #     cbar.set_label('Density contrast ($g/cm^{3}$)', fontsize=fs, labelpad=-1)
    # ******************************************************************************************************************

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ------------------------------------------------------------------------------------------------------------------
    # SAVE FIGURE TO DISC
    if use_tight_layout:
        if file_type == "svg":
            plt.savefig(file_path+'.'+file_type, bbox_inches='tight', dpi=720, format='svg')
        elif file_type == "pdf":
            plt.savefig(file_path+'.'+file_type, bbox_inches='tight', dpi=720, format='pdf')
        elif file_type == "ps":
            plt.savefig(file_path+'.'+file_type, bbox_inches='tight', dpi=720, format='ps')
        elif file_type == "eps":
            plt.savefig(file_path+'.'+file_type, bbox_inches='tight', dpi=720, format='eps')
        elif file_type == "png":
            plt.savefig(file_path+'.'+file_type, bbox_inches='tight', dpi=720, format='png')
    else:
        if file_type == "svg":
            plt.savefig(file_path+'.'+file_type, dpi=720, format='svg')
        elif file_type == "pdf":
            plt.savefig(file_path+'.'+file_type, dpi=720, format='pdf')
        elif file_type == "ps":
            plt.savefig(file_path+'.'+file_type, dpi=720, format='ps')
        elif file_type == "eps":
            plt.savefig(file_path+'.'+file_type, dpi=720, format='eps')
        elif file_type == "png":
            plt.savefig(file_path+'.'+file_type, dpi=720, format='png')
    # ------------------------------------------------------------------------------------------------------------------
    return
