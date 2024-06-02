"""
SAVE GMG MODEL AS A EDITABLE VECTOR GRAPHICS FIGURE.

CALLED FROM gmg draw_model() FUNC
"""

#import pygmt
import numpy as np, copy
import matplotlib as mpl
import os

def plot_fig(file_path, file_type, use_tight_layout, fs, 
             aspect_ratio, ps, calc_line_width, font_type,
             t_frame_min, t_frame_max, 
             grav_frame_min, grav_frame_max, 
             mag_frame_min, mag_frame_max,
             draw_polygons, polygon_alpha, 
             draw_fixed_layers, layer_line_width,
             draw_floating_layers, layer_line_alpha,
             draw_colorbar, colorbar_x, colorbar_y, 
             colorbar_size_x, colorbar_size_y, 
             draw_xy_data, xy_size, xy_color,
             draw_wells, well_fs, well_line_width, 
             draw_faults, faults_lw, layer_list, fault_list,
             observed_topography_list, observed_gravity_list, 
             observed_magnetic_list, 
             outcrop_data_list, well_data_list,
             segy_data_list, t_frame, g_frame, m_frame,
             predicted_gravity, gravity_rms_value,
             predicted_nt, magnetic_rms_value, 
             area, xp):
    """
    CREATE A FIGURE FROM THE CURRENT GMG MODEL USING USER SPECIFIED 
    PLOTTING OPTIONS AND SAVE TO AN EXTERNAL FILE
    """

    # -------------------------------------------------------------------------
    # SET DEFAULT PLOTTING PARAMS

    # FIGURE FONT SIZE
    fs = fs
    
    # MARKER SIZE
    ps = ps

    # LINE WIDTH
    lw = calc_line_width

    y_shift_value = "10c"

    # DIR CONTAINING BOREHOLE ICON
    print(os.path.dirname(os.path.abspath(__file__)))
    borehole_dir = os.path.dirname(os.path.abspath(__file__))+'../docs/icons/'
    print(borehole_dir)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # START FIGURE PLOTTING

    fig = pygmt.Figure()

    # SET THE MODEL AXES
    x_start_model, x_end_model, z_end_model, z_start_model = np.array(area)/1000

   
    # ---------------------------------------------------------------------------
    # *****************************************************************************
    # PLOT MODEL PANEL
    # *****************************************************************************

    # PLOT BASEMAP
    fig.basemap(region=[x_start_model, x_end_model, z_end_model, z_start_model], 
                    projection="x1/"+str(aspect_ratio),
                    frame=["SWne", "xaf+lDistance (km)", "yaf+lDepth (km)"])
    
    # # PLOT SEGY
    # for segy in segy_data_list:
    #     if segy.mpl_actor.get_visible() is True:
    #         fig.

    # PLOT LAYER POLYGONS
    if draw_polygons is True:
        print(layer_list)
        for layer in layer_list[1:]:

            # GET ACTOR ATTRIBUTES
            x = layer.polygon_mpl_actor[0].get_xy()[:,0]
            y = layer.polygon_mpl_actor[0].get_xy()[:,1]
            fc =  layer.polygon_mpl_actor[0].get_fc()
            # CONVERT POLYGON COLOR TO RGB
            fc_rgb = np.array([fc[0]*255, fc[1]*255, fc[2]*255]).astype(int)
            fa =  layer.polygon_mpl_actor[0].get_alpha()*100
            lc = layer.color
            print(x)
            print(y)
            print(fc_rgb)
            print(fa)
            print(lc)
            print(layer_line_width)

            # DRAW LAYER FILL
            fig.plot(x=x, y=y, transparency=fa, fill=str(str(fc_rgb[0])+"/"+str(fc_rgb[1])+"/"+str(fc_rgb[2])))

            # DRAW LAYER LINE
            fig.plot(x=x, y=y, pen=str(layer_line_width)+"p,", transparency=fa)
            
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

    fig.shift_origin(yshift=y_shift_value)
     # ---------------------------------------------------------------------------
    # # PLOT TOPOGRAPHY PANEL
    if t_frame.get_visible() is True:

        fig.basemap(region=[x_start_model, x_end_model, t_frame_min, t_frame_max], 
                    projection="x1/"+str(aspect_ratio),
                    frame=["SWne", "xaf+lDistance (km)", "yaf+lMag anom (nT)"])

    # ---------------------------------------------------------------------------
    
    # fig.shift_origin(yshift=y_shift_value)

    # ---------------------------------------------------------------------------
    # # PLOT GRAVITY PANEL
    # if g_frame.get_visible() is True:
    # ---------------------------------------------------------------------------

    # fig.shift_origin(yshift=y_shift_value)

    # ---------------------------------------------------------------------------
    # # PLOT MAGNETIC PANEL
    # if m_frame.get_visible() is True:
    #     ax3 = plt.subplot2grid((num_rows, 1), (row_counter, 0), 
    #                               rowspan=1, colspan=1)
    #     row_counter += 1

    #     # AXIS OPTIONS
    #     plt.ylabel('nT', fontsize=fs)
    #     plt.xlim(x_start_model, x_end_model)
    #     plt.ylim(mag_frame_min, mag_frame_max)
    #     ax3.set_xticks([])

    #     # PLOT OBSERVED MAGNETIC DATA
    #     if observed_magnetic_list is not None:
    #         for i in range(0, len(observed_magnetic_list)):
    #             if observed_magnetic_list[i] is not None:
    #                 ax2.scatter(observed_magnetic_list[i].data[:, 0], 
    #                   observed_magnetic_list[i].data[:, 0], marker='o',
    #                     color=observed_magnetic_list[i].color, s=ps)

    #     # PLOT PREDICTED MAGNETIC ANOMALY
    #     if predicted_nt is not None:
    #         ax2.plot(xp * 0.001, predicted_nt, color='g', 
    #                  linewidth=calc_line_width)

    #     # PLOT RMS VALUE
    #     if magnetic_rms_value is not None:
    #         ax2.annotate('RMS misfit: ' + str(magnetic_rms_value), xy=(10, 60), 
    #            xytext=(10, 60), fontsize=fs, horizontalalignment='left', 
    #               clip_on=False)

    #     # SET AXIS POSITIONS
    #     ax3.spines['right'].set_color('none')
    #     ax3.spines['bottom'].set_position('center')
    #     ax3.spines['top'].set_position('center')
    #     ax3.tick_params(axis='x', which='both', labelbottom='off', 
    #                     labeltop='off')
    #     ax3.tick_params(axis='y', which='both', left='on', 
    #                       right='off', labelright='off')
    # ---------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # SAVE FIGURE TO DISC
        fig.savefig(file_path+'.'+file_type, dpi=720)
        # ------------------------------------------------------------------------------------------------------------------
    return
