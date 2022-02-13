from __future__ import print_function
import os
import sys
import argparse
import plot_util as p
from colors import *


parser = argparse.ArgumentParser(description="plot settings")
parser.add_argument('-i' , '--input_dir'              , dest='input_dir'       , help='input dir path '                        , required=True      )
parser.add_argument('-o' , '--output_dir'             , dest='output_dir'      , help='output dir path'                        , default='plots/'                                                 )
parser.add_argument('-n' , '--nbins'                  , dest='nbins'           , help='number of bins for the histograms'      , default=30                                                      )
parser.add_argument('-xt', '--xaxis_title'            , dest='xaxis_title'     , help='X-axis title'                           , default=None                                                    )
parser.add_argument('-yt', '--yaxis_title'            , dest='yaxis_title'     , help='X-axis title'                           , default=None                                                    )
parser.add_argument('-y' , '--yaxis_range'            , dest='yaxis_range'     , help='Y-axis range set by user'               , default=None                                                    )
parser.add_argument('-x' , '--xaxis_range'            , dest='xaxis_range'     , help='X-axis range set by user'               , default=None                                                    )
parser.add_argument('-l' , '--yaxis_log'              , dest='yaxis_log'       , help='Y-axis set to log'                      , default=False                              , action='store_true')
parser.add_argument('-s' , '--sig_scale'              , dest='sig_scale'       , help='Signal scale'                           , default=1                                                       )
parser.add_argument('-b' , '--bkg_scale'              , dest='bkg_scale'       , help='Background scale'                       , default=1                                                       )
parser.add_argument('-u' , '--rm_udflow'              , dest='rm_udflow'       , help='Remove underflow'                       , default=False                              , action='store_true')
parser.add_argument('-ov' , '--rm_ovflow'              , dest='rm_ovflow'       , help='Remove overflow'                       , default=False                              , action='store_true')
parser.add_argument('-au' , '--add_udflow'              , dest='add_udflow'       , help='Add underflow'                       , default=False                              , action='store_true')
parser.add_argument('-aov' , '--add_ovflow'              , dest='add_ovflow'       , help='Add overflow'                       , default=True                              , action='store_true')
parser.add_argument('-d' , '--draw_data'              , dest='draw_data'       , help='Draw data'                              , default=False                              , action='store_true')
parser.add_argument('-1' , '--stack_signal'           , dest='stack_signal'    , help='stack signal'                           , default=True                              , action='store_true')
parser.add_argument('-O' , '--output_name'            , dest='output_name'     , help='output file name when plot single hist' , default=None                                                    )
parser.add_argument('-w' , '--whatSR'                 , dest='whatSR'          , help='what selecton for the nine bins'        , default="FatJetsSDMassCut"                                      )
parser.add_argument('-R' , '--right_hand'             , dest='right_hand'      , help='remove right side (<)'                , default=False                      , action='store_true')
parser.add_argument('-OP' , '--do_optimize'           , dest='do_optimize'     , help='do cut optimization'                , default=False                      , action='store_true')
parser.add_argument('-F' , '--do_fit'           , dest='do_fit'     , help='do fit'                , default=False                      , action='store_true')

parser.add_argument('hist_name', metavar='<histogram_names>=(e.g. FatJetsSDMassCut__hh_pt)', type=str, nargs='*', help='patterns to use to filter histograms to dump')

args = parser.parse_args()


hist_name = None
if len(args.hist_name) == 0:
    print("")
    print("---------------------------------------------------------------")
    print("No hist name provided, I will plot all histograms in region "+args.whatSR)
    print("If you want to plot something else, please provide a histogram name.")
    print(" e.g. FatJetsSDMassCut__hh_pt")
    print("---------------------------------------------------------------")
    print("")
else:
    hist_name = ','.join(args.hist_name)


sig_scale = float(args.sig_scale)
bkg_scale = float(args.bkg_scale)
right_hand = args.right_hand
output_name = args.output_name
do_optimize = args.do_optimize
do_fit = args.do_fit

input_dir = args.input_dir
if input_dir[-1] != "/":
    input_dir += "/"

output_dir = args.output_dir + "/" + input_dir.split("/")[-3] + "/" +input_dir.split("/")[-2]
if do_optimize:
    output_dir = args.output_dir + "/" + input_dir.split("/")[-3]+"_cutOptimize" + "/" + input_dir.split("/")[-2]

os.system("mkdir -p "+args.output_dir)
os.system("cp index.php "+args.output_dir)
os.system("mkdir -p "+args.output_dir + "/" + input_dir.split("/")[-3] )
os.system("cp index.php "+args.output_dir + "/" + input_dir.split("/")[-3] )
os.system("mkdir -p "+output_dir)
os.system("cp index.php "+output_dir)

year = -1
lumi = -1
if "2016" in input_dir:
    year = 2016
    lumi = 35.9
elif "2017" in input_dir:
    year = 2017
    lumi = 41.3
elif "2018" in input_dir:
    year = 2018
    lumi = 59.7
else:
    year = 6051
    lumi = 137

#sig_fnames = [input_dir+"ttbar.root"]
sig_fnames = [input_dir+"tt1L.root"]
#bkg_fnames = [input_dir+"others.root", input_dir+"Higgs.root", input_dir+"VH.root", input_dir+"ttH.root", input_dir+"tW.root", input_dir+"qcd.root", input_dir+"ttbar.root"]
#bkg_fnames = [input_dir+"others.root", input_dir+"Higgs.root", input_dir+"VH.root", input_dir+"ttH.root", input_dir+"qcd.root", input_dir+"ttbar.root"]
#bkg_fnames = [input_dir+"others.root", input_dir+"Higgs.root", input_dir+"VH.root", input_dir+"ttH.root", input_dir+"ttbar.root", input_dir+"qcd.root"]
#bkg_fnames = [input_dir+"WJets.root", input_dir+"qcd.root", input_dir+"tt1L.root", input_dir+"tt2L.root"]
bkg_fnames = [input_dir+"others.root", input_dir+"qcd.root", input_dir+"tt2L.root"]

#sig_legends = ["t#bar{t}+jets"]
#sig_legends = ["t#bar{t}+jets"]
sig_legends = ["t#bar{t}+jets 1L"]
#bkg_legends = ["others", "ggH+VBFH", "VH", "t#bar{t}H", "tW", "QCD", "t#bar{t}+jets"]
#bkg_legends = ["others", "ggH+VBFH", "VH", "t#bar{t}H", "QCD", "t#bar{t}+jets"]
#bkg_legends = ["others", "ggH+VBFH", "VH", "t#bar{t}H", "t#bar{t}+jets", "QCD"]
#bkg_legends = ["w(lv)+jets", "V+jets,VV", "QCD", "t#bar{t}+jets 1L", "t#bar{t}+jets 2L"]
bkg_legends = ["V+jets,VV", "QCD", "t#bar{t}+jets 2L"]


data_fname = input_dir+"data.root"
sig_colors = [2005,617, 839, 800, 1, 632]
#bkg_colors = [2001, 2003, 2011, 920, 2007, 46, 2005, 800]
#bkg_colors = [2001, 2003, 2011, 920, 2007, 2005, 800, 839]
#bkg_colors = [2001, 2003, 2011, 920, 2005, 2007, 800, 839]
bkg_colors = [2001, 2007, 800, 839]


if float(args.sig_scale) != 1:
    for idx in range(len(sig_legends)):
        sig_legends[idx] = sig_legends[idx] + " x "+args.sig_scale


if hist_name:
    if do_optimize:
        print("optimize cut for "+hist_name)
        p.makeplot_cutOptimize(
            sig_fnames_=sig_fnames,
            bkg_fnames_=bkg_fnames,
            sig_legends_=sig_legends,
            bkg_legends_=bkg_legends,
            sig_colors_=sig_colors,
            bkg_colors_=bkg_colors,
            hist_name_=hist_name,
            sig_scale_=sig_scale,
            bkg_scale_=bkg_scale,
            dir_name_=output_dir,
            right_hand_=right_hand,
            output_name_=args.output_name,
            extraoptions={
                "remove_underflow":args.rm_udflow,
                "remove_overflow":args.rm_ovflow,
                "lumi_value": lumi,
                "xaxis_label": args.xaxis_title,
                },
            )

    elif do_fit:
        print("fit for "+hist_name)
        p.makeplot_fitmass(
            data_fname_=data_fname,
            hist_name_=hist_name,
            dir_name_=output_dir,
            output_name_=args.output_name,
            extraoptions={
                "nbins":int(args.nbins),
                "remove_underflow":args.rm_udflow,
                "remove_overflow":args.rm_ovflow,
                "lumi_value": lumi,
                "xaxis_label": args.xaxis_title,
                },
            )
    else:
        print("plotting "+hist_name)
        p.makeplot_single(
            sig_fnames_=sig_fnames,
            bkg_fnames_=bkg_fnames,
            data_fname_=data_fname if args.draw_data else None,
            sig_legends_=sig_legends,
            bkg_legends_=bkg_legends,
            sig_colors_=sig_colors,
            bkg_colors_=bkg_colors,
            hist_name_=hist_name, 
            sig_scale_=float(args.sig_scale),
            dir_name_=output_dir, 
            output_name_=args.output_name,
            extraoptions={ 
                "nbins":int(args.nbins),
                "yaxis_range":args.yaxis_range.split(',') if args.yaxis_range else [],
                "xaxis_range":args.xaxis_range.split(',') if args.xaxis_range else [],
                "xaxis_label": args.xaxis_title,
                "yaxis_label": args.yaxis_title,
                "remove_underflow":args.rm_udflow,
                "remove_overflow":args.rm_ovflow,
                "add_underflow":args.add_udflow,
                "add_overflow":args.add_ovflow,
                "lumi_value": lumi,
                "ratio_range": [0., 2.],
                "xaxis_label": args.xaxis_title,
                "stack_signal": args.stack_signal
                },
            )
            
else:
    print("plotting everything in "+args.whatSR)
    SR_all= [args.whatSR]
    p.makeplot_all(
        sig_fnames_=sig_fnames,
        bkg_fnames_=bkg_fnames,
        data_fname_=data_fname if args.draw_data else None,
        sig_legends_=sig_legends,
        bkg_legends_=bkg_legends,
        sig_colors_=sig_colors,
        bkg_colors_=bkg_colors,
        regions_name_=SR_all, 
        sig_scale_=float(args.sig_scale),
        dir_name_=output_dir, 
        extraoptions={ 
            "nbins":int(args.nbins),
            "yaxis_range":args.yaxis_range.split(',') if args.yaxis_range else [],
            "xaxis_range":args.xaxis_range.split(',') if args.xaxis_range else [],
            "xaxis_label": args.xaxis_title,
            "yaxis_label": args.yaxis_title,
            "remove_underflow":args.rm_udflow,
            "remove_overflow":args.rm_ovflow,
            "add_underflow":args.add_udflow,
            "add_overflow":args.add_ovflow,
            "lumi_value": lumi,
            "ratio_range": [0.6, 1.2],
            "stack_signal": args.stack_signal
            },
        )
 



