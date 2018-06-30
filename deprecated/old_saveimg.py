def save_key_img(dv,hk,efv,key,modelmap=None,weightmap=None,
                 component='Bulk',filename="Map.fits"):

    myhdr=dv[key].maps.header
    weights = dv[key].maps.wts
    if weightmap != None:
        weights=weightmap
    else:
        if hk.cfp.mask == True:
            bi=np.where(dv[key].mapping.arcmap > hk.cfp.maskrad)
            weights[bi] = 0.0
            if hk.hk_ins.name == 'abell_2146':
                bi = np.where(dv[key].mapping.angmap < hk.hk_ins.lban_nw)
                weights[bi] = 0.0
                bi = np.where(dv[key].mapping.angmap > hk.hk_ins.uban_nw)
                weights[bi] = 0.0

    tstr = 'Test_Run_'
    if hk.cfp.testmode == False:
        tstr = 'Full_Run_'
        
    if modelmap==None:

        if component == 'Bulk' or component == 'Residual':
            modelsky = mlf.get_best_map(efv,hk,dv)
            beam_conv = ip.conv_inst_beam(modelsky,dv[key].mapping.pixsize,instrument=dv[key].mapping.instrument)
            gc_model=ip.apply_xfer(beam_conv, dv[key].mapping.tab,instrument=dv[key].mapping.instrument)
            model = gc_model

        if component == 'Bulk':
    ### Let's try using the "Header Date Unit" (HDU) format/style.
            hdu1 = fits.PrimaryHDU(modelsky,header=myhdr)
#            hdu1.header = dv[key].mapping.w.to_header()
            hdu1.header.append(("Title","Jansky/beam Sky Map"))
            hdu1.header.append(("Target",hk.hk_ins.name))
            hdu2 = fits.ImageHDU(beam_conv)
            hdu2.header = dv[key].mapping.w.to_header()
            hdu2.name = 'BeamConvolved'
            hdu2.header.append(("Title","Beam Convolved Map"))
            hdu2.header.append(("Target",hk.hk_ins.name))
            hdu2.header.append(("XTENSION","What Mate"))
            hdu2.header.append(("SIMPLE","T")) # DAMN, THIS SHIT REAL
            hdu2.verify('fix')
            hdu3 = fits.ImageHDU(gc_model)
            hdu3.header = dv[key].mapping.w.to_header()
            hdu3.header.append(("Title","Fully Filtered Map"))
            hdu3.header.append(("Target",hk.hk_ins.name))
            hdu3.header.append(("XTENSION","What Mate"))
            hdu3.name = 'FilteredMap';
            hdu3.header.append(("SIMPLE","T")) # DAMN, THIS SHIT REAL
            hdu3.verify('fix')
            hdu4 = fits.ImageHDU(weights)
            hdu4.header = dv[key].mapping.w.to_header()
            hdu4.header.append(("Title","Weight Map"))
            hdu4.header.append(("Target",hk.hk_ins.name))
            hdu4.header.append(("XTENSION","What Mate"))
            hdu4.header.append(("SIMPLE","T")) # DAMN, THIS SHIT REAL
            hdu4.name = 'WeightMap'
            hdu4.verify('fix')
            for i in range(len(efv.values)):
                hdu1.header.append(("Bin_"+str(i),efv.values[i],"keV cm**-3"))
            
            hdu1.header.add_history("Fits were performed on "+hk.log.today+ ".")
            hdulist = fits.HDUList([hdu1,hdu2,hdu3,hdu4])
            hdulist.info()
            filename=tstr+"Fitted_bulk_cluster.fits"
            fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
            hdulist.writeto(fullpath,overwrite=True,output_verify="exception")
#            import pdb; pdb.set_trace()
#        fits.writeto(fullpath,ptsrc_map,ptsrchdr,overwrite=True)
#        fits.update(fullpath,dv[key].maps.wts,1)

        if component == 'Ptsrc' or component == 'Residual':
            ptsrc_map, ptsrc_hdr = fits.getdata(hk.cfp.psfile, header=True)
            if hk.cfp.ndim > hk.cfp.bins:
                addind=0
                if hk.cfp.mn_lvl == True:
                    addind+=1
                if hk.cfp.ptsrc == True:
                    myfac = efv.psolns[hk.cfp.bins+addind].value
                    addind+=1
                else:
                    raise Exception("No point source was fit dummy!")
#        import pdb; pdb.set_trace()
            ptsrc_map *= myfac[0]
            ps_model=ip.apply_xfer(ptsrc_map, dv[key].mapping.tab,instrument=dv[key].mapping.instrument)
            model += ps_model
                
        if component == 'Residual':
            residual = dv[key].maps.data - model
            hdu1 = fits.PrimaryHDU(residual,header=myhdr)
#            hdu1.header = dv[key].mapping.w.to_header()
            hdu1.header.append(("Title","Residual Map"))
            hdu1.header.append(("Target",hk.hk_ins.name))
            hdu2 = fits.ImageHDU(dv[key].maps.data)
            hdu2.header = dv[key].mapping.w.to_header()
            hdu2.header.append(("Title","Data Map"))
            hdu2.header.append(("Target",hk.hk_ins.name))
            hdu2.verify('fix')
            hdu3 = fits.ImageHDU(model)
            hdu3.header = dv[key].mapping.w.to_header()
            hdu3.header.append(("Title","Model (all components included)"))
            hdu3.header.append(("Target",hk.hk_ins.name))
            hdu3.verify('fix')
            hdulist = fits.HDUList([hdu1,hdu2,hdu3])        
            filename=tstr+"Residual.fits"
            fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
            hdulist.writeto(fullpath,overwrite=True)
 #       fits.writeto(fullpath,ptsrc_map,ptsrchdr,overwrite=True)
 #       fits.update(fullpath,dv[key].maps.wts,1)

### Now, the loop for if modelmap == something...: 
    else:
        hdu1 = fits.PrimaryHDU(modelmap,header=myhdr)
#        hdu1.header = dv[key].mapping.w.to_header()
        hdu1.header.append(("Title","Jansky/beam Sky Map"))
        hdu1.header.append(("Target",hk.hk_ins.name))
        hdulist = fits.HDUList([hdu1])
        if weightmap != None:
            hdu2 = fits.ImageHDU(weightmap)
            hdu2.header = dv[key].mapping.w.to_header()
            hdu2.header.append(("Title","Weight Map"))
            hdu2.header.append(("Target",hk.hk_ins.name))
            hdu2.verify('fix')
            hdulist = fits.HDUList([hdu1,hdu2])
            
        fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+tstr+filename)
#        print fullpath,hdu1.header
        hdu1.header.append(("Simple","T"))
        hdu1.verify('fix')
#        import pdb; pdb.set_trace()
        hdulist.writeto(fullpath,overwrite=True)
