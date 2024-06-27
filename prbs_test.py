#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Prbs Test
# GNU Radio version: v3.8.5.0-6-g57bd109d

from distutils.version import StrictVersion

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print("Warning: failed to XInitThreads()")

from PyQt5 import Qt
from gnuradio import eng_notation
from gnuradio import qtgui
from gnuradio.filter import firdes
import sip
from cmath import exp, pi; from random import randint; import numpy as np;
from gnuradio import channels
from gnuradio import digital
from gnuradio import filter
from gnuradio import gr
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio.qtgui import Range, RangeWidget
import mapper

from gnuradio import qtgui

class prbs_test(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Prbs Test")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Prbs Test")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "prbs_test")

        try:
            if StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
                self.restoreGeometry(self.settings.value("geometry").toByteArray())
            else:
                self.restoreGeometry(self.settings.value("geometry"))
        except:
            pass

        ##################################################
        # Variables
        ##################################################
        self.sro = sro = 0
        self.samp_rate = samp_rate = 2000000
        self.samp_per_sym = samp_per_sym = 2
        self.pream = pream = (mapper.preamble_generator(32,511,1033)).get_preamble()
        self.frame_width = frame_width = 1024*2
        self.cfo = cfo = 0
        self.SNR = SNR = 12

        ##################################################
        # Blocks
        ##################################################
        self._sro_range = Range(-1e-1, 1e-1, 1e-3, 0, 200)
        self._sro_win = RangeWidget(self._sro_range, self.set_sro, 'Sample Rate Offset', "counter_slider", float)
        self.top_layout.addWidget(self._sro_win)
        self._cfo_range = Range(-1e-1, 1e-1, 1e-3, 0, 200)
        self._cfo_win = RangeWidget(self._cfo_range, self.set_cfo, 'Freq Offset', "counter_slider", float)
        self.top_layout.addWidget(self._cfo_win)
        self._SNR_tool_bar = Qt.QToolBar(self)
        self._SNR_tool_bar.addWidget(Qt.QLabel('SNR' + ": "))
        self._SNR_line_edit = Qt.QLineEdit(str(self.SNR))
        self._SNR_tool_bar.addWidget(self._SNR_line_edit)
        self._SNR_line_edit.returnPressed.connect(
            lambda: self.set_SNR(eng_notation.str_to_num(str(self._SNR_line_edit.text()))))
        self.top_layout.addWidget(self._SNR_tool_bar)
        self.root_raised_cosine_filter_0_0 = filter.interp_fir_filter_ccf(
            samp_per_sym*0+1,
            firdes.root_raised_cosine(
                1,
                samp_rate,
                samp_rate/samp_per_sym,
                0.35,
                11))
        self.root_raised_cosine_filter_0 = filter.interp_fir_filter_ccf(
            samp_per_sym,
            firdes.root_raised_cosine(
                1,
                samp_rate,
                samp_rate/samp_per_sym,
                0.35,
                11))
        self.qtgui_time_sink_x_0 = qtgui.time_sink_c(
            1024, #size
            samp_rate, #samp_rate
            "", #name
            1 #number of inputs
        )
        self.qtgui_time_sink_x_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0.set_y_axis(-1, 1)

        self.qtgui_time_sink_x_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0.enable_tags(True)
        self.qtgui_time_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0.enable_autoscale(True)
        self.qtgui_time_sink_x_0.enable_grid(False)
        self.qtgui_time_sink_x_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0.enable_stem_plot(False)


        labels = ['Signal 1', 'Signal 2', 'Signal 3', 'Signal 4', 'Signal 5',
            'Signal 6', 'Signal 7', 'Signal 8', 'Signal 9', 'Signal 10']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ['blue', 'red', 'green', 'black', 'cyan',
            'magenta', 'yellow', 'dark red', 'dark green', 'dark blue']
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]
        styles = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1]


        for i in range(2):
            if len(labels[i]) == 0:
                if (i % 2 == 0):
                    self.qtgui_time_sink_x_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_win)
        self.qtgui_const_sink_x_0_0 = qtgui.const_sink_c(
            840-84, #size
            'QT GUI Plot', #name
            1 #number of inputs
        )
        self.qtgui_const_sink_x_0_0.set_update_time(0.10)
        self.qtgui_const_sink_x_0_0.set_y_axis(-0.0015, 0.0015)
        self.qtgui_const_sink_x_0_0.set_x_axis(-0.0015, 0.0015)
        self.qtgui_const_sink_x_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, "")
        self.qtgui_const_sink_x_0_0.enable_autoscale(False)
        self.qtgui_const_sink_x_0_0.enable_grid(False)
        self.qtgui_const_sink_x_0_0.enable_axis_labels(True)


        labels = ['', '', '', '', '',
            '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ["blue", "red", "red", "red", "red",
            "red", "red", "red", "red", "red"]
        styles = [0, 0, 0, 0, 0,
            0, 0, 0, 0, 0]
        markers = [0, 0, 0, 0, 0,
            0, 0, 0, 0, 0]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_const_sink_x_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_const_sink_x_0_0.set_line_label(i, labels[i])
            self.qtgui_const_sink_x_0_0.set_line_width(i, widths[i])
            self.qtgui_const_sink_x_0_0.set_line_color(i, colors[i])
            self.qtgui_const_sink_x_0_0.set_line_style(i, styles[i])
            self.qtgui_const_sink_x_0_0.set_line_marker(i, markers[i])
            self.qtgui_const_sink_x_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_const_sink_x_0_0_win = sip.wrapinstance(self.qtgui_const_sink_x_0_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_const_sink_x_0_0_win)
        self.mapper_preamble_sync_demapper_0 = mapper.preamble_sync_demapper(frame_width, pream, mapper.BPSK, [0,1], 1, 3, False)
        self.mapper_preamble_insert_bb_0 = mapper.preamble_insert_bb(frame_width, pream)
        self.mapper_prbs_source_b_0 = mapper.prbs_source_b("PRBS31", frame_width-len(pream))
        self.mapper_prbs_sink_b_0 = mapper.prbs_sink_b("PRBS31", frame_width-len(pream))
        self.mapper_mapper_0 = mapper.mapper(mapper.BPSK, [0,1])
        self.digital_costas_loop_cc_0 = digital.costas_loop_cc(1e-2, 2, False)
        self.digital_clock_recovery_mm_xx_0 = digital.clock_recovery_mm_cc(samp_per_sym*(1+0.0), 0.25*0.175*0.175, 0.5, 0.175, 0.005)
        self.digital_binary_slicer_fb_0 = digital.binary_slicer_fb()
        self.channels_channel_model_0 = channels.channel_model(
            noise_voltage=(10**(-SNR/10.0))*np.sqrt(samp_per_sym),
            frequency_offset=cfo,
            epsilon=1.0+sro,
            taps=[1],
            noise_seed=0,
            block_tags=False)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.channels_channel_model_0, 0), (self.root_raised_cosine_filter_0_0, 0))
        self.connect((self.digital_binary_slicer_fb_0, 0), (self.mapper_prbs_sink_b_0, 0))
        self.connect((self.digital_clock_recovery_mm_xx_0, 0), (self.mapper_preamble_sync_demapper_0, 0))
        self.connect((self.digital_clock_recovery_mm_xx_0, 0), (self.qtgui_const_sink_x_0_0, 0))
        self.connect((self.digital_costas_loop_cc_0, 0), (self.digital_clock_recovery_mm_xx_0, 0))
        self.connect((self.mapper_mapper_0, 0), (self.root_raised_cosine_filter_0, 0))
        self.connect((self.mapper_prbs_source_b_0, 0), (self.mapper_preamble_insert_bb_0, 0))
        self.connect((self.mapper_preamble_insert_bb_0, 0), (self.mapper_mapper_0, 0))
        self.connect((self.mapper_preamble_sync_demapper_0, 0), (self.digital_binary_slicer_fb_0, 0))
        self.connect((self.root_raised_cosine_filter_0, 0), (self.channels_channel_model_0, 0))
        self.connect((self.root_raised_cosine_filter_0_0, 0), (self.digital_costas_loop_cc_0, 0))
        self.connect((self.root_raised_cosine_filter_0_0, 0), (self.qtgui_time_sink_x_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "prbs_test")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_sro(self):
        return self.sro

    def set_sro(self, sro):
        self.sro = sro
        self.channels_channel_model_0.set_timing_offset(1.0+self.sro)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.qtgui_time_sink_x_0.set_samp_rate(self.samp_rate)
        self.root_raised_cosine_filter_0.set_taps(firdes.root_raised_cosine(1, self.samp_rate, self.samp_rate/self.samp_per_sym, 0.35, 11))
        self.root_raised_cosine_filter_0_0.set_taps(firdes.root_raised_cosine(1, self.samp_rate, self.samp_rate/self.samp_per_sym, 0.35, 11))

    def get_samp_per_sym(self):
        return self.samp_per_sym

    def set_samp_per_sym(self, samp_per_sym):
        self.samp_per_sym = samp_per_sym
        self.channels_channel_model_0.set_noise_voltage((10**(-self.SNR/10.0))*np.sqrt(self.samp_per_sym))
        self.digital_clock_recovery_mm_xx_0.set_omega(self.samp_per_sym*(1+0.0))
        self.root_raised_cosine_filter_0.set_taps(firdes.root_raised_cosine(1, self.samp_rate, self.samp_rate/self.samp_per_sym, 0.35, 11))
        self.root_raised_cosine_filter_0_0.set_taps(firdes.root_raised_cosine(1, self.samp_rate, self.samp_rate/self.samp_per_sym, 0.35, 11))

    def get_pream(self):
        return self.pream

    def set_pream(self, pream):
        self.pream = pream
        self.mapper_preamble_insert_bb_0.mapper.set_preamble(self.pream)

    def get_frame_width(self):
        return self.frame_width

    def set_frame_width(self, frame_width):
        self.frame_width = frame_width
        self.mapper_preamble_insert_bb_0.mapper.set_width(self.frame_width)

    def get_cfo(self):
        return self.cfo

    def set_cfo(self, cfo):
        self.cfo = cfo
        self.channels_channel_model_0.set_frequency_offset(self.cfo)

    def get_SNR(self):
        return self.SNR

    def set_SNR(self, SNR):
        self.SNR = SNR
        Qt.QMetaObject.invokeMethod(self._SNR_line_edit, "setText", Qt.Q_ARG("QString", eng_notation.num_to_str(self.SNR)))
        self.channels_channel_model_0.set_noise_voltage((10**(-self.SNR/10.0))*np.sqrt(self.samp_per_sym))





def main(top_block_cls=prbs_test, options=None):

    if StrictVersion("4.5.0") <= StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    def quitting():
        tb.stop()
        tb.wait()

    qapp.aboutToQuit.connect(quitting)
    qapp.exec_()

if __name__ == '__main__':
    main()
