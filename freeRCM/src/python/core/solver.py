import ctypes
import os
import numpy as np
import thermodynamics
import data_structures

orgProps = thermodynamics.orgProps
dict2struct = data_structures.dict2struct


class Params(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.POINTER(ctypes.c_double)),
        ("P", ctypes.c_double),
        ("antProps", ctypes.POINTER(ctypes.c_double)),
        ("nrtlA", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("nrtlB", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("nrtlC", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
        ("TcCel", ctypes.POINTER(ctypes.c_double)),
        ("Pc", ctypes.POINTER(ctypes.c_double)),
        ("omega", ctypes.POINTER(ctypes.c_double)),
        ("Ncomps", ctypes.c_int),
        ("antMethod", ctypes.c_int),
        ("actMethod", ctypes.c_int),
        ("dxi", ctypes.c_double),
        ("n_it", ctypes.c_int),
        ("maxiter", ctypes.c_int),
        ("ftol", ctypes.c_double),
        ("xtol", ctypes.c_double),
    ]


class Curves(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.POINTER(ctypes.c_double)),
        ("y", ctypes.POINTER(ctypes.c_double)),
        ("T", ctypes.POINTER(ctypes.c_double)),
    ]


def RCM(comps, selected_comps, P, allProps, opts, x0n, genOpt):

    try:
        dxi = float(opts.dxi) * 1.0
    except TypeError:
        print("type err\n")
        dxi = 0.02
    except AttributeError:
        print("attr err\n")
        dxi = 0.02

    try:
        n_it = opts.n_it * 1
    except TypeError:
        n_it = 250
    except AttributeError:
        n_it = 250

    Ncomps = len(selected_comps)
    if opts.antMethod == 1:
        props = orgProps(1, comps, selected_comps, allProps)
    else:
        props = orgProps(3, comps, selected_comps, allProps)

    if opts.activity == 2:
        props2 = orgProps(2, comps, selected_comps, allProps)
        for key in props2:
            props[key] = props2[key]

    if opts.activity == 3:
        props2 = orgProps(2, comps, selected_comps, allProps)
        props3 = orgProps(4, comps, selected_comps, allProps)
        for key in props2:
            props[key] = props2[key]
        for key in props3:
            props[key] = props3[key]

    # Load native libraries with proper dependency resolution
    lib_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "lib")

    # Set library search path for macOS/Linux
    if os.name == 'posix':
        current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
        os.environ['LD_LIBRARY_PATH'] = f"{lib_dir}:{current_ld_path}"
    elif os.name == 'nt':
        current_path = os.environ.get('PATH', '')
        os.environ['PATH'] = f"{lib_dir};{current_path}"

    # Load dependencies first, then main library
    minpack_path = os.path.join(lib_dir, "libminpack.so")
    solver_path = os.path.join(lib_dir, "RCM_solver.so")

    # Load in dependency order
    ctypes.CDLL(minpack_path)
    lib = ctypes.CDLL(solver_path)

    if genOpt == 2:
        Nlines = 1
        x0 = x0n
    else:
        Nlines = opts.lines
        x_bank = np.linspace(0.27, 0.49, opts.lines)

    xPlot = np.zeros([2 * n_it, Ncomps, Nlines])
    yPlot = np.zeros([2 * n_it, Ncomps, Nlines])
    TPlot = np.zeros([2 * n_it, 1, Nlines])

    Ncomps_c = (ctypes.c_int)(Ncomps)
    P_c = (ctypes.c_double)(P)
    TcCel_c = np.asfortranarray(props.TcCel).ctypes.data_as(
        ctypes.POINTER(ctypes.c_double)
    )
    Pc_c = np.asfortranarray(props.Pc).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    omega_c = np.asfortranarray(props.omega).ctypes.data_as(
        ctypes.POINTER(ctypes.c_double)
    )
    nrtlA_c = np.asfortranarray(props.NRTL_aij).ctypes.data_as(
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
    )
    nrtlB_c = np.asfortranarray(props.NRTL_bij).ctypes.data_as(
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
    )
    nrtlC_c = np.asfortranarray(props.NRTL_cij).ctypes.data_as(
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
    )

    if opts.antMethod == 1:
        antProps_c = props.antoine.flatten(order="C").ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)
        )
    else:
        antProps_c = props.PLXANT.flatten(order="C").ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)
        )
    dxi_c = (ctypes.c_double)(dxi)
    n_it_c = (ctypes.c_int)(n_it)
    maxiter_c = (ctypes.c_int)(opts.lmopts["maxiter"])
    ftol_c = (ctypes.c_double)(opts.lmopts["ftol"])
    xtol_c = (ctypes.c_double)(opts.lmopts["xtol"])

    params = Params(
        P=P_c,
        antProps=antProps_c,
        nrtlA=nrtlA_c,
        nrtlB=nrtlB_c,
        nrtlC=nrtlC_c,
        TcCel=TcCel_c,
        Pc=Pc_c,
        omega=omega_c,
        Ncomps=Ncomps_c,
        antMethod=opts.antMethod,
        actMethod=opts.activity,
        dxi=dxi_c,
        n_it=n_it_c,
        maxiter=maxiter_c,
        ftol=ftol_c,
        xtol=xtol_c,
    )
    lib.RCM.argtypes = [ctypes.POINTER(Params)]
    lib.RCM.restype = ctypes.POINTER(Curves)

    for k in range(0, Nlines):

        if genOpt != 2:
            x0 = np.array([x_bank[k], 1.5 - 3 * x_bank[k], -0.5 + 2 * x_bank[k]])

        x_c = x0.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        params.x = x_c

        curves = lib.RCM(ctypes.byref(params))

        xPlot[:, :, k] = np.ctypeslib.as_array(
            curves.contents.x, shape=(2 * n_it, Ncomps)
        )
        yPlot[:, :, k] = np.ctypeslib.as_array(
            curves.contents.y, shape=(2 * n_it, Ncomps)
        )
        TPlot[:, :, k] = np.ctypeslib.as_array(curves.contents.T, shape=(2 * n_it, 1))

        lib.freeCurveMem(curves)
        del curves

    allCurves = dict2struct()
    allCurves.x = xPlot
    allCurves.y = yPlot
    allCurves.T = TPlot

    return allCurves
