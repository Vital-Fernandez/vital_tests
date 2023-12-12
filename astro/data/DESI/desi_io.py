


def read_spectra(
        infile,
        single=False,
        targetids=None,
        rows=None,
        skip_hdus=None,
        select_columns={
            "FIBERMAP": None,
            "EXP_FIBERMAP": None,
            "SCORES": None,
            "EXTRA_CATALOG": None,
        },
):
    """
    Read Spectra object from FITS file.

    This reads data written by the write_spectra function.  A new Spectra
    object is instantiated and returned.

    Args:
        infile (str): path to read
        single (bool): if True, keep spectra as single precision in memory.
        targetids (list): Optional, list of targetids to read from file, if present.
        rows (list): Optional, list of rows to read from file
        skip_hdus (list): Optional, list/set/tuple of HDUs to skip
        select_columns (dict): Optional, dictionary to select column names to be read. Default, all columns are read.

    Returns (Spectra):
        The object containing the data read from disk.

    `skip_hdus` options are FIBERMAP, EXP_FIBERMAP, SCORES, EXTRA_CATALOG, MASK, RESOLUTION;
    where MASK and RESOLUTION mean to skip those for all cameras.
    Note that WAVE, FLUX, and IVAR are always required.

    If a table HDU is not listed in `select_columns`, all of its columns will be read

    User can optionally specify targetids OR rows, but not both
    """
    log = get_logger()
    infile = checkgzip(infile)
    ftype = np.float64
    if single:
        ftype = np.float32

    infile = os.path.abspath(infile)
    if not os.path.isfile(infile):
        raise IOError("{} is not a file".format(infile))

    t0 = time.time()
    hdus = fitsio.FITS(infile, mode="r")
    nhdu = len(hdus)

    if targetids is not None and rows is not None:
        raise ValueError('Set rows or targetids but not both')

    # - default skip_hdus empty set -> include everything, without
    # - having to check for None before checking if X is in skip_hdus
    if skip_hdus is None:
        skip_hdus = set()

    # - Map targets -> rows and exp_rows.
    # - Note: coadds can have uncoadded EXP_FIBERMAP HDU with more rows than
    # - the coadded FIBERMAP HDU, so track rows vs. exp_rows separately
    exp_rows = None
    if targetids is not None:
        targetids = np.atleast_1d(targetids)
        file_targetids = hdus["FIBERMAP"].read(columns="TARGETID")
        rows = np.where(np.isin(file_targetids, targetids))[0]
        if 'EXP_FIBERMAP' in hdus and 'EXP_FIBERMAP' not in skip_hdus:
            exp_targetids = hdus["EXP_FIBERMAP"].read(columns="TARGETID")
            exp_rows = np.where(np.isin(exp_targetids, targetids))[0]
        if len(rows) == 0:
            return Spectra()
    elif rows is not None:
        rows = np.asarray(rows)
        # figure out exp_rows
        file_targetids = hdus["FIBERMAP"].read(rows=rows, columns="TARGETID")
        if 'EXP_FIBERMAP' in hdus and 'EXP_FIBERMAP' not in skip_hdus:
            exp_targetids = hdus["EXP_FIBERMAP"].read(columns="TARGETID")
            exp_rows = np.where(np.isin(exp_targetids, file_targetids))[0]

    if select_columns is None:
        select_columns = dict()

    for extname in ("FIBERMAP", "EXP_FIBERMAP", "SCORES", "EXTRA_CATALOG"):
        if extname not in select_columns:
            select_columns[extname] = None

    # load the metadata.
    meta = dict(hdus[0].read_header())

    # initialize data objects

    bands = []
    fmap = None
    expfmap = None
    wave = None
    flux = None
    ivar = None
    mask = None
    res = None
    extra = None
    extra_catalog = None
    scores = None

    # For efficiency, go through the HDUs in disk-order.  Use the
    # extension name to determine where to put the data.  We don't
    # explicitly copy the data, since that will be done when constructing
    # the Spectra object.

    for h in range(1, nhdu):
        name = hdus[h].read_header()["EXTNAME"]
        log.debug('Reading %s', name)
        if name == "FIBERMAP":
            if name not in skip_hdus:
                fmap = encode_table(
                    Table(
                        hdus[h].read(rows=rows, columns=select_columns["FIBERMAP"]),
                        copy=True,
                    ).as_array()
                )
        elif name == "EXP_FIBERMAP":
            if name not in skip_hdus:
                expfmap = encode_table(
                    Table(
                        hdus[h].read(rows=exp_rows, columns=select_columns["EXP_FIBERMAP"]),
                        copy=True,
                    ).as_array()
                )
        elif name == "SCORES":
            if name not in skip_hdus:
                scores = encode_table(
                    Table(
                        hdus[h].read(rows=rows, columns=select_columns["SCORES"]),
                        copy=True,
                    ).as_array()
                )
        elif name == "EXTRA_CATALOG":
            if name not in skip_hdus:
                extra_catalog = encode_table(
                    Table(
                        hdus[h].read(
                            rows=rows, columns=select_columns["EXTRA_CATALOG"]
                        ),
                        copy=True,
                    ).as_array()
                )
        else:
            # Find the band based on the name
            mat = re.match(r"(.*)_(.*)", name)
            if mat is None:
                raise RuntimeError(
                    "FITS extension name {} does not contain the band".format(name)
                )
            band = mat.group(1).lower()
            type = mat.group(2)
            if band not in bands:
                bands.append(band)
            if type == "WAVELENGTH":
                if wave is None:
                    wave = {}
                # - Note: keep original float64 resolution for wavelength
                wave[band] = native_endian(hdus[h].read())
            elif type == "FLUX":
                if flux is None:
                    flux = {}
                flux[band] = _read_image(hdus, h, ftype, rows=rows)
            elif type == "IVAR":
                if ivar is None:
                    ivar = {}
                ivar[band] = _read_image(hdus, h, ftype, rows=rows)
            elif type == "MASK" and type not in skip_hdus:
                if mask is None:
                    mask = {}
                mask[band] = _read_image(hdus, h, np.uint32, rows=rows)
            elif type == "RESOLUTION" and type not in skip_hdus:
                if res is None:
                    res = {}
                res[band] = _read_image(hdus, h, ftype, rows=rows)
            elif type != "MASK" and type != "RESOLUTION" and type not in skip_hdus:
                # this must be an "extra" HDU
                log.debug('Reading extra HDU %s', name)
                if extra is None:
                    extra = {}
                if band not in extra:
                    extra[band] = {}

                extra[band][type] = _read_image(hdus, h, ftype, rows=rows)

    hdus.close()
    duration = time.time() - t0
    log.info(iotime.format("read", infile, duration))

    # Construct the Spectra object from the data.  If there are any
    # inconsistencies in the sizes of the arrays read from the file,
    # they will be caught by the constructor.

    spec = Spectra(
        bands,
        wave,
        flux,
        ivar,
        mask=mask,
        resolution_data=res,
        fibermap=fmap,
        exp_fibermap=expfmap,
        meta=meta,
        extra=extra,
        extra_catalog=extra_catalog,
        single=single,
        scores=scores,
    )

    # if needed, sort spectra to match order of targetids, which could be
    # different than the order they appear in the file
    if targetids is not None:
        from desispec.util import ordered_unique
        # - Input targetids that we found in the file, in the order they appear in targetids
        ii = np.isin(targetids, spec.fibermap['TARGETID'])
        found_targetids = ordered_unique(targetids[ii])
        log.debug('found_targetids=%s', found_targetids)

        # - Unique targetids of input file in the order they first appear
        input_targetids = ordered_unique(spec.fibermap['TARGETID'])
        log.debug('input_targetids=%s', np.asarray(input_targetids))

        # - Only reorder if needed
        if not np.all(input_targetids == found_targetids):
            rows = np.concatenate([np.where(spec.fibermap['TARGETID'] == tid)[0] for tid in targetids])
            log.debug("spec.fibermap['TARGETID'] = %s", np.asarray(spec.fibermap['TARGETID']))
            log.debug("rows for subselection=%s", rows)
            spec = spec[rows]

    return spec