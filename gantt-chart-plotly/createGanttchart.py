# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 11:52:42 2018

@author: vlachos
"""

import plotly.plotly as py
import plotly.figure_factory as ff
import plotly

df = [# Phase 1 - Orientation
      dict(Task='Orientation and Proposal', Start='2018-10-08', Finish='2018-11-04', Resource='Complete'),
      dict(Task='Brief Literature Review', Start='2018-10-08', Finish='2018-10-26', Resource='Complete'),
      dict(Task='Research Proposal', Start='2018-10-22', Finish='2018-11-04', Resource='Complete'),
      # Phase 2 - Data Collection
      dict(Task='Data Collection/Algorithms', Start='2018-11-05', Finish='2018-12-02', Resource='Incomplete'),
      # Data Preparation
      dict(Task='Data Preparation algorithms', Start='2018-11-05', Finish='2018-12-09', Resource='Incomplete'),
      # Phase 3 - Exploration/Preprocessing
      dict(Task='Exploration-Preprocessing', Start='2018-12-03', Finish='2019-01-13', Resource='Incomplete'),

      # Phase 4a - Data Analysis
      dict(Task='Machine Learning', Start='2018-12-24', Finish='2019-03-31', Resource='Not Started'),
      # Important Milestones
      dict(Task='Kick-Off meeting', Start='2018-12-03', Finish='2018-12-04', Resource='Milestone'),
      dict(Task='Send minutes', Start='2018-12-05', Finish='2018-12-06', Resource='Milestone'),
      dict(Task='Progress Meeting 1', Start='2019-01-14', Finish='2019-01-15', Resource='Milestone'),
      dict(Task='Submit Report (contents etc.)', Start='2019-02-01', Finish='2019-02-02', Resource='Milestone'),
      dict(Task='Mid-Term Review', Start='2019-02-18', Finish='2019-02-19', Resource='Milestone'),
      dict(Task='Send minutes', Start='2019-02-20', Finish='2019-02-21', Resource='Milestone'),
      dict(Task='Progress Meeting 2', Start='2019-03-18', Finish='2019-03-19', Resource='Milestone'),
      dict(Task='Submit Draft Report', Start='2019-05-06', Finish='2019-05-07', Resource='Milestone'),
      dict(Task='Submit Final Report', Start='2019-05-20', Finish='2019-05-21', Resource='Milestone'),
      # Report Writing
      dict(Task='Introduction', Start='2018-12-03', Finish='2019-01-31', Resource='Report Writing'),
      dict(Task='Methodology', Start='2019-01-07', Finish='2019-04-07', Resource='Report Writing'),
      dict(Task='Results', Start='2019-01-21', Finish='2019-04-14', Resource='Report Writing'),
      dict(Task='Discussion', Start='2019-02-25', Finish='2019-05-05', Resource='Report Writing')
      ]

colors = {'Not Started': 'rgb(191, 191, 191)',
          'Incomplete': 'rgb(204, 204, 0)',
          'Complete': 'rgb(0, 255, 100)',
          'Milestone': 'rgb(255, 0, 0)',
          'Report Writing': 'rgb(153, 153, 255)'}


fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True)
fig['layout'].update(autosize=False, width=3508/2, height=2480/2, margin=dict(l=200))
plotly.offline.plot(fig, filename='gantt-group-tasks-together')#, world_readable=True)

