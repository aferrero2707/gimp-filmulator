#include <stdlib.h>
#include <gegl.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <libgimpbase/gimpbase.h>
#include <libgimpwidgets/gimpwidgets.h>
#include <glib/gi18n.h>
#include <string.h>

#include <iostream>
#include <fstream>

#include <lcms2.h>
#include <gexiv2/gexiv2-metadata.h>

#include <tiffio.h>

#define VERSION "0.0.1"

//#define HAVE_GIMP_2_9 1


static cmsCIExyYTRIPLE srgb_primaries = {
    {0.6400, 0.3300, 1.0},
    {0.3000, 0.6000, 1.0},
    {0.1500, 0.0600, 1.0}
};

static cmsCIExyYTRIPLE srgb_primaries_pre_quantized = {
    {0.639998686, 0.330010138, 1.0},
    {0.300003784, 0.600003357, 1.0},
    {0.150002046, 0.059997204, 1.0}
};

static cmsCIExyYTRIPLE rec2020_primaries = {
    {0.7079, 0.2920, 1.0},
    {0.1702, 0.7965, 1.0},
    {0.1314, 0.0459, 1.0}
};

static cmsCIExyYTRIPLE rec2020_primaries_prequantized = {
    {0.708012540607, 0.291993664388, 1.0},
    {0.169991652439, 0.797007778423, 1.0},
    {0.130997824007, 0.045996550894, 1.0}
};

/* ************************** WHITE POINTS ************************** */

/* D65 WHITE POINTS */

static cmsCIExyY  d65_srgb_adobe_specs = {0.3127, 0.3290, 1.0};
/* White point from the sRGB.icm and AdobeRGB1998 profile specs:
 * http://www.adobe.com/digitalimag/pdfs/AdobeRGB1998.pdf
 * 4.2.1 Reference Display White Point
 * The chromaticity coordinates of white displayed on
 * the reference color monitor shall be x=0.3127, y=0.3290.
 * . . . [which] correspond to CIE Standard Illuminant D65.
 *
 * Wikipedia gives this same white point for SMPTE-C.
 * This white point is also given in the sRGB color space specs.
 * It's probably correct for most or all of the standard D65 profiles.
 *
 * The D65 white point values used in the LCMS virtual sRGB profile
 * is slightly different than the D65 white point values given in the
 * sRGB color space specs, so the LCMS virtual sRGB profile
 * doesn't match sRGB profiles made using the values given in the
 * sRGB color space specs.
 *
 * */

static cmsHPROFILE make_srgb_profile()
{
  /* ***** Make profile: sRGB, D65, sRGB TRC */
  /*
   * */
  cmsCIExyYTRIPLE primaries = srgb_primaries_pre_quantized;
  cmsCIExyY whitepoint = d65_srgb_adobe_specs;
  cmsToneCurve* tone_curve[3];
  cmsHPROFILE profile = NULL;

  /* sRGB TRC */
  cmsFloat64Number srgb_parameters[5] =
  { 2.4, 1.0 / 1.055,  0.055 / 1.055, 1.0 / 12.92, 0.04045 };
  cmsToneCurve *curve = cmsBuildParametricToneCurve(NULL, 4, srgb_parameters);
  tone_curve[0] = tone_curve[1] = tone_curve[2] = curve;

  profile = cmsCreateRGBProfile ( &whitepoint, &primaries, tone_curve );
  cmsMLU *copyright = cmsMLUalloc(NULL, 1);
  cmsMLUsetASCII(copyright, "en", "US", "Copyright 2015, Elle Stone (website: http://ninedegreesbelow.com/; email: ellestone@ninedegreesbelow.com). This ICC profile is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License (https://creativecommons.org/licenses/by-sa/3.0/legalcode).");
  cmsWriteTag(profile, cmsSigCopyrightTag, copyright);
  // V4
  cmsMLU *description = cmsMLUalloc(NULL, 1);
  cmsMLUsetASCII(description, "en", "US", "sRGB-elle-V4.icc");
  cmsWriteTag(profile, cmsSigProfileDescriptionTag, description);
  //const char* filename = "Rec2020-elle-V4-rec709.icc";
  //cmsSaveProfileToFile(profile, filename);
  cmsMLUfree(description);
  /**/

  return( profile );
}



static cmsHPROFILE make_rec2020_linear_profile()
{
  /* ***** Make profile: Rec.2020, D65, Rec709 TRC */
  /*
   * */
  cmsCIExyYTRIPLE primaries = rec2020_primaries_prequantized;
  cmsCIExyY whitepoint = d65_srgb_adobe_specs;
  /* rec.709 */
  cmsToneCurve* tone_curve[3];
  cmsHPROFILE profile = NULL;

  cmsToneCurve *curve = cmsBuildGamma (NULL, 1.00);
  tone_curve[0] = tone_curve[1] = tone_curve[2] = curve;

  profile = cmsCreateRGBProfile ( &whitepoint, &primaries, tone_curve );
  cmsMLU *copyright = cmsMLUalloc(NULL, 1);
  cmsMLUsetASCII(copyright, "en", "US", "Copyright 2015, Elle Stone (website: http://ninedegreesbelow.com/; email: ellestone@ninedegreesbelow.com). This ICC profile is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License (https://creativecommons.org/licenses/by-sa/3.0/legalcode).");
  cmsWriteTag(profile, cmsSigCopyrightTag, copyright);
  // V4
  cmsMLU *description = cmsMLUalloc(NULL, 1);
  cmsMLUsetASCII(description, "en", "US", "Rec2020-elle-V4.icc");
  cmsWriteTag(profile, cmsSigProfileDescriptionTag, description);
  //const char* filename = "Rec2020-elle-V4-rec709.icc";
  //cmsSaveProfileToFile(profile, filename);
  cmsMLUfree(description);
  /**/

  return( profile );
}



static cmsHPROFILE make_rec2020_srgb_profile()
{
  /* ***** Make profile: Rec.2020, D65, Rec709 TRC */
  /*
   * */
  cmsCIExyYTRIPLE primaries = rec2020_primaries_prequantized;
  cmsCIExyY whitepoint = d65_srgb_adobe_specs;
  /* rec.709 */
  cmsToneCurve* tone_curve[3];
  cmsHPROFILE profile = NULL;

  /* sRGB TRC */
  cmsFloat64Number srgb_parameters[5] =
  { 2.4, 1.0 / 1.055,  0.055 / 1.055, 1.0 / 12.92, 0.04045 };
  cmsToneCurve *curve = cmsBuildParametricToneCurve(NULL, 4, srgb_parameters);
  tone_curve[0] = tone_curve[1] = tone_curve[2] = curve;

  profile = cmsCreateRGBProfile ( &whitepoint, &primaries, tone_curve );
  cmsMLU *copyright = cmsMLUalloc(NULL, 1);
  cmsMLUsetASCII(copyright, "en", "US", "Copyright 2015, Elle Stone (website: http://ninedegreesbelow.com/; email: ellestone@ninedegreesbelow.com). This ICC profile is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License (https://creativecommons.org/licenses/by-sa/3.0/legalcode).");
  cmsWriteTag(profile, cmsSigCopyrightTag, copyright);
  // V4
  cmsMLU *description = cmsMLUalloc(NULL, 1);
  cmsMLUsetASCII(description, "en", "US", "Rec2020-elle-V4.icc");
  cmsWriteTag(profile, cmsSigProfileDescriptionTag, description);
  //const char* filename = "Rec2020-elle-V4-rec709.icc";
  //cmsSaveProfileToFile(profile, filename);
  cmsMLUfree(description);
  /**/

  return( profile );
}



// Manage different versions of the GIMP API.
#define _gimp_item_is_valid gimp_item_is_valid
#define _gimp_image_get_item_position gimp_image_get_item_position
#if GIMP_MINOR_VERSION<=8
#define _gimp_item_get_visible gimp_drawable_get_visible
#else
#define _gimp_item_get_visible gimp_item_get_visible
#endif

void query();
void run(const gchar *name,
    gint nparams,
    const GimpParam *param,
    gint *nreturn_vals,
    GimpParam **return_vals);

//long pf_save_gimp_image(ufraw_data *uf, GtkWidget *widget);

GimpPlugInInfo PLUG_IN_INFO = {
    NULL,  /* init_procedure */
    NULL,  /* quit_procedure */
    query, /* query_procedure */
    run,   /* run_procedure */
};

MAIN()

void query()
{
  static const GimpParamDef args[] = {
      {GIMP_PDB_INT32,    (gchar*)"run_mode", (gchar*)"Interactive, non-interactive"},
      {GIMP_PDB_IMAGE,    (gchar*)"image",    (gchar*)"Input image"},
      {GIMP_PDB_DRAWABLE, (gchar*)"drawable", (gchar*)"Input drawable (unused)"},
  };

  gimp_install_procedure("plug-in-filmulator",             // name
      "Filmulator",                    // blurb
      "Filmulator",                    // help
      "Carlo Vaccari", // author
      "Carlo Vaccari", // copyright
      "2016",                     // date
      "_Filmulator...",                // menu_path
      "RGB*",              // image_types
      GIMP_PLUGIN,                // type
      G_N_ELEMENTS(args),         // nparams
      0,                          // nreturn_vals
      args,                       // params
      0);                         // return_vals

  gimp_plugin_menu_register("plug-in-filmulator", "<Image>/Filters");
}

char *pf_binary;
gboolean sendToGimpMode;


void run(const gchar *name,
    gint nparams,
    const GimpParam *param,
    gint *nreturn_vals,
    GimpParam **return_vals)
{
  GimpRunMode run_mode = (GimpRunMode)param[0].data.d_int32;

  int size;
  GimpPDBStatusType status = GIMP_PDB_SUCCESS;

#if !GLIB_CHECK_VERSION(2,31,0)
  g_thread_init(NULL);
#endif
#ifndef _WIN32
  gdk_threads_init();
  gdk_threads_enter();
#endif

#if HAVE_GIMP_2_9
  gegl_init(NULL, NULL);
#endif

#if HAVE_GIMP_2_9
  GeglBuffer *buffer;
#else
  GimpDrawable *drawable;
  GimpPixelRgn pixel_region;
  int tile_height, row, nrows;
#endif
  //gint32 layer;

  static GimpParam return_values[1];
  *return_vals  = return_values;
  *nreturn_vals = 1;
  return_values[0].type = GIMP_PDB_STATUS;

  int image_id = param[1].data.d_drawable;
#if GIMP_MINOR_VERSION<=8
  gimp_tile_cache_ntiles(2*(gimp_image_width(image_id)/gimp_tile_width() + 1));
#endif

  // Retrieve the list of desired layers.
  int nb_layers = 0;
  int *layers = gimp_image_get_layers(image_id,&nb_layers),
      active_layer_id = gimp_image_get_active_layer(image_id);
  int source_layer_id = active_layer_id;

  GimpParasite *filmulator_parasite = gimp_item_get_parasite( active_layer_id, "filmulator-config" );
  std::cout<<"Filmulator plug-in: filmulator_parasite="<<filmulator_parasite<<std::endl;
  bool replace_layer = false;
  std::string pfiname;
  if( filmulator_parasite && gimp_parasite_data_size( filmulator_parasite ) > 0 &&
      gimp_parasite_data( filmulator_parasite ) != NULL ) {
    replace_layer = true;
  }

  if( replace_layer ) {
    if (active_layer_id>=0) {
      int i = 0; for (i = 0; i<nb_layers; ++i) if (layers[i]==active_layer_id) break;
      if (i<nb_layers - 1) source_layer_id = layers[i + 1];
      else source_layer_id = -1;
    }
  }
  std::cout<<"Filmulator plug-in: replace_layer="<<replace_layer<<"  source_layer_id="<<source_layer_id<<std::endl;

  //GimpParasite *exif_parasite = gimp_image_parasite_find( image_id, "gimp-image-metadata" );
  GimpMetadata* exif_metadata = gimp_image_get_metadata( image_id );
  GimpParasite *icc_parasite  = gimp_image_parasite_find( image_id, "icc-profile" );
  glong iccsize = 0;
  void* iccdata = NULL;

  std::cout<<std::endl<<std::endl
      <<"image_id: "<<image_id
      <<"  ICC parasite: "<<icc_parasite
      <<"  EXIF metadata: "<<exif_metadata
      <<std::endl<<std::endl;

  cmsHPROFILE gimp_profile = NULL;
  cmsHTRANSFORM transform1 = NULL;
  cmsHTRANSFORM transform2 = NULL;

  if( icc_parasite && gimp_parasite_data_size( icc_parasite ) > 0 &&
      gimp_parasite_data( icc_parasite ) != NULL ) {
    iccsize = gimp_parasite_data_size( icc_parasite );
    iccdata = malloc( iccsize );
    memcpy( iccdata, gimp_parasite_data( icc_parasite ), iccsize );
    gimp_profile = cmsOpenProfileFromMem( iccdata, iccsize );
  } else {
    gimp_profile = make_srgb_profile();
  }

  cmsHPROFILE rec2020_lin = make_rec2020_linear_profile();
  cmsHPROFILE rec2020_srgb = make_rec2020_srgb_profile();

  if( gimp_profile && rec2020_lin && rec2020_srgb ) {
    transform1 = cmsCreateTransform( gimp_profile, TYPE_RGB_FLT,
        rec2020_lin, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC,
        cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );

    transform2 = cmsCreateTransform( rec2020_srgb, TYPE_RGB_FLT,
        gimp_profile, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC,
        cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );
  }

  if(gimp_profile) cmsCloseProfile( gimp_profile );
  if(rec2020_lin) cmsCloseProfile( rec2020_lin );
  if(rec2020_srgb) cmsCloseProfile( rec2020_srgb );

  std::string filename;
  cmsBool is_lin_gamma = false;
  std::string format = "R'G'B' float";

  if( source_layer_id >= 0 ) {
    // Get input buffer
    gint rgn_x, rgn_y, rgn_width, rgn_height;
    if (!_gimp_item_is_valid(source_layer_id)) return;
    if (!gimp_drawable_mask_intersect(source_layer_id,&rgn_x,&rgn_y,&rgn_width,&rgn_height)) return;
    const int spectrum = (gimp_drawable_is_rgb(source_layer_id)?3:1) +
        (gimp_drawable_has_alpha(source_layer_id)?1:0);
    unsigned char* inbuf = (unsigned char*)malloc( sizeof(float)*3*rgn_width*rgn_height );
#if GIMP_MINOR_VERSION<=8
    GimpDrawable *drawable = gimp_drawable_get(source_layer_id);
    GimpPixelRgn region;
    gimp_pixel_rgn_init(&region,drawable,rgn_x,rgn_y,rgn_width,rgn_height,false,false);
    guchar* row = g_new(guchar,rgn_width*spectrum), *ptrs = 0;
    switch (spectrum) {
    case 1 : {
      float* ptr_r = (float*)inbuf;
      for( int y = 0; y < rgn_height; y++ ) {
        gimp_pixel_rgn_get_row(&region,ptrs=row,rgn_x,rgn_y + y,rgn_width);
        for( int x = 0; x < rgn_width; x++ ) {
          float val = ((float)*(ptrs++))/255;
          *(ptr_r++) = val;
          *(ptr_r++) = val;
          *(ptr_r++) = val;
        }
      }
    } break;
    case 2 : {
      float* ptr_r = (float*)inbuf;
      for( int y = 0; y < rgn_height; y++ ) {
        gimp_pixel_rgn_get_row(&region,ptrs=row,rgn_x,rgn_y + y,rgn_width);
        for( int x = 0; x < rgn_width; x++ ) {
          float val = ((float)*(ptrs++))/255; ptrs++;
          *(ptr_r++) = val;
          *(ptr_r++) = val;
          *(ptr_r++) = val;
        }
      }
    } break;
    case 3 : {
      float* ptr_r = (float*)inbuf;
      for( int y = 0; y < rgn_height; y++ ) {
        gimp_pixel_rgn_get_row(&region,ptrs=row,rgn_x,rgn_y + y,rgn_width);
        for( int x = 0; x < rgn_width; x++ ) {
          *(ptr_r++) = ((float)*(ptrs++))/255;
          *(ptr_r++) = ((float)*(ptrs++))/255;
          *(ptr_r++) = ((float)*(ptrs++))/255;
        }
      }
    } break;
    case 4 : {
      float* ptr_r = (float*)inbuf;
      for( int y = 0; y < rgn_height; y++ ) {
        gimp_pixel_rgn_get_row(&region,ptrs=row,rgn_x,rgn_y + y,rgn_width);
        for( int x = 0; x < rgn_width; x++ ) {
          *(ptr_r++) = ((float)*(ptrs++))/255;
          *(ptr_r++) = ((float)*(ptrs++))/255;
          *(ptr_r++) = ((float)*(ptrs++))/255; ptrs++;
        }
      }
    } break;
    }
    g_free(row);
    gimp_drawable_detach(drawable);
#else
    if( iccdata ) {
      cmsHPROFILE iccprofile = cmsOpenProfileFromMem( iccdata, iccsize );
      if( iccprofile ) {
        char tstr[1024];
        cmsGetProfileInfoASCII(iccprofile, cmsInfoDescription, "en", "US", tstr, 1024);
        std::cout<<std::endl<<std::endl<<"embedded profile: "<<tstr<<std::endl<<std::endl;
        cmsToneCurve *red_trc   = (cmsToneCurve*)cmsReadTag(iccprofile, cmsSigRedTRCTag);
        is_lin_gamma = cmsIsToneCurveLinear(red_trc);
        cmsCloseProfile( iccprofile );
      }
    }

    GeglRectangle rect;
    gegl_rectangle_set(&rect,rgn_x,rgn_y,rgn_width,rgn_height);
    buffer = gimp_drawable_get_buffer(source_layer_id);
    format = is_lin_gamma ? "RGB float" : "R'G'B' float";
    gegl_buffer_get(buffer,&rect,1,babl_format(format.c_str()),inbuf,0,GEGL_ABYSS_NONE);
    g_object_unref(buffer);
#endif

    if( transform1 )
      cmsDoTransform( transform1, inbuf, inbuf, rgn_width*rgn_height );

    // Write layer buffer into TIFF file on disk
    std::string filename = "/tmp/gimp-filmulator-in.tif";

    int sampleperpixel = 3;

    TIFF *out= TIFFOpen(filename.c_str(), "w");
    TIFFSetField (out, TIFFTAG_IMAGEWIDTH, rgn_width);  // set the width of the image
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, rgn_height);    // set the height of the image
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, sizeof(float)*8);    // set the size of the channels
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
    //   Some other essential fields to set that you do not have to understand for now.
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    tsize_t linebytes = sampleperpixel * rgn_width * sizeof(float);     // length in memory of one row of pixel in the image.

    std::cout<<"Writing output TIFF file..."<<std::endl;
    std::cout<<"  linebytes="<<linebytes<<"  TIFFScanlineSize(out)="<<TIFFScanlineSize(out)<<std::endl;
    unsigned char *tiffbuf = NULL;        // buffer used to store the row of pixel information for writing to file
    //    Allocating memory to store the pixels of current row
    if (TIFFScanlineSize(out) < linebytes)
      tiffbuf =(unsigned char *)_TIFFmalloc(linebytes);
    else
      tiffbuf = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(out));

    // We set the strip size of the file to be size of one row of pixels
    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, sizeof(float)*rgn_width*sampleperpixel));

    std::cout<<"  inbuf="<<(void*)inbuf<<"  tiffbuf="<<(void*)tiffbuf<<std::endl;
    //Now writing image to the file one strip at a time
    for (uint32 row = 0; row < rgn_height; row++) {
      //std::cout<<"  writing row "<<row<<"..."<<std::endl;
      memcpy(tiffbuf, &inbuf[(rgn_height-row-1)*linebytes], linebytes);
      //std::cout<<"  buffer copied"<<std::endl;
      if (TIFFWriteScanline(out, tiffbuf, row, 0) < 0)
        break;
      //std::cout<<"  ... done"<<std::endl;
    }
    std::cout<<"... done"<<std::endl;

    (void) TIFFClose(out);

    if (tiffbuf)
      _TIFFfree(tiffbuf);

    free( inbuf );

    if( exif_metadata ) {
      gexiv2_metadata_save_file( /*(GExiv2Metadata*)*/exif_metadata, filename.c_str(), NULL );
      g_object_unref( exif_metadata );
    }

    std::cout<<"plug-in: run_mode="<<run_mode<<"  GIMP_RUN_INTERACTIVE="<<GIMP_RUN_INTERACTIVE<<std::endl;

    /* BUG - what should be done with GIMP_RUN_WITH_LAST_VALS */
    if (run_mode == GIMP_RUN_INTERACTIVE) {

      // run filmulator and pass the temporray tiff file to it
      std::string outfilename = "/tmp/gimp-filmulator-out.tif";

      char filmulator_cmd[1000];
      //snprintf( filmulator_cmd, 999, "filmulator \"%s\" \"%s\"", filename.c_str(), outfilename.c_str() );
      snprintf( filmulator_cmd, 999, "cp \"%s\" \"%s\"", filename.c_str(), outfilename.c_str() );
      std::cout<<"Running filmulator: "<<filmulator_cmd<<std::endl;
      system( filmulator_cmd );
      std::cout<<"Filmulator finished"<<std::endl;

      unlink( filename.c_str() );

      uint32 width, height;

      TIFF *tiffin= TIFFOpen(outfilename.c_str(), "r");
      TIFFGetField(tiffin, TIFFTAG_IMAGEWIDTH, &width);           // uint32 width;
      TIFFGetField(tiffin, TIFFTAG_IMAGELENGTH, &height);        // uint32 height;
      uint32 npixels=width*height;
      unsigned char* raster=(unsigned char *) _TIFFmalloc(npixels *sizeof(float)*sampleperpixel);

      if( raster ) {
        // Transfer the output layers back into GIMP.
        std::cout<<"PhF plug-in: copying buffer..."<<std::endl;
        GimpLayerModeEffects layer_blendmode = GIMP_NORMAL_MODE;
        gint layer_posx = 0, layer_posy = 0;
        double layer_opacity = 100;

        gint32 dest_layer_id = active_layer_id;
        if( !replace_layer ) {
          /* Create a new layer to hold the output... */
          gint32 layer = gimp_layer_new(image_id, _("PhF output"), width,
              height, GIMP_RGB_IMAGE, 100.0,
              GIMP_NORMAL_MODE);
          std::cout<<"PhF plug-in: new layer created"<<std::endl;
#if defined(GIMP_CHECK_VERSION) && GIMP_CHECK_VERSION(2,7,3)
          gimp_image_insert_layer(image_id, layer, 0, -1);
#else
          gimp_image_add_layer(image_id, layer, -1);
#endif
          std::cout<<"Filmulator plug-in: new layer added"<<std::endl;
          dest_layer_id = layer;
        }
        /* Get the drawable and set the pixel region for our load... */
#if HAVE_GIMP_2_9
        buffer = gimp_drawable_get_buffer(dest_layer_id);
#else
        drawable = gimp_drawable_get(dest_layer_id);
        gimp_pixel_rgn_init(&pixel_region, drawable, 0, 0, drawable->width,
            drawable->height, TRUE, FALSE);
        tile_height = gimp_tile_height();
#endif

        //TIFFReadRGBImage(tiffin, width, height, raster, 0);

        linebytes = sampleperpixel * width * sizeof(float);     // length in memory of one row of pixel in the image.

        tiffbuf = NULL;        // buffer used to store the row of pixel information for writing to file
        //    Allocating memory to store the pixels of current row
        if (TIFFScanlineSize(tiffin) < linebytes)
          tiffbuf =(unsigned char *)_TIFFmalloc(linebytes);
        else
          tiffbuf = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(tiffin));

        //Now reading image from the file one strip at a time
        for (uint32 row = 0; row < height; row++) {
          if (TIFFReadScanline(tiffin, tiffbuf, row, 0) < 0)
            break;
          memcpy(&raster[(height-row-1)*linebytes], tiffbuf, linebytes);
        }

        (void) TIFFClose(tiffin);

        if (tiffbuf)
          _TIFFfree(tiffbuf);


        if( transform2 )
          cmsDoTransform( transform2, raster, raster, npixels );

        unlink( outfilename.c_str() );

#if HAVE_GIMP_2_9
        format = is_lin_gamma ? "RGB float" : "R'G'B' float";
        GeglRectangle gegl_rect;
        gegl_rect.x = 0;
        gegl_rect.y = 0;
        gegl_rect.width = width;
        gegl_rect.height = height;
        gegl_buffer_set(buffer, &gegl_rect,
            //GEGL_RECTANGLE(0, 0, width, height),
            0, babl_format(format.c_str()), raster,
            GEGL_AUTO_ROWSTRIDE);
        //gimp_drawable_merge_shadow(layer_id,true);
        gimp_drawable_update(dest_layer_id,0,0,width,height);
        gimp_layer_resize(dest_layer_id,width,height,0,0);
#else
        for (row = 0; row < Crop.height; row += tile_height) {
          nrows = MIN(Crop.height - row, tile_height);
          gimp_pixel_rgn_set_rect(&pixel_region,
              uf->thumb.buffer + 3 * row * Crop.width, 0, row, Crop.width, nrows);
        }
#endif
        std::cout<<"Filmulator plug-in: buffer copied"<<std::endl;

        _TIFFfree( raster );

        //GimpParasite *cfg_parasite;
        //cfg_parasite = gimp_parasite_new("filmulator-config",
        //    GIMP_PARASITE_PERSISTENT, strlen(buffer), buffer);
        //gimp_item_attach_parasite(dest_layer_id, cfg_parasite);
        //gimp_parasite_free(cfg_parasite);
      }
    }

#if HAVE_GIMP_2_9
    gegl_buffer_flush(buffer);
    g_object_unref(buffer);
#else
    gimp_drawable_flush(drawable);
    gimp_drawable_detach(drawable);
#endif

  }

  if( transform1 ) cmsDeleteTransform( transform1 );
  if( transform2 ) cmsDeleteTransform( transform2 );

  gimp_displays_flush();

  std::cout<<"Plug-in: setting return values"<<std::endl;
  return_values[0].data.d_status = status;
  std::cout<<"Plug-in: return values done"<<std::endl;
  std::cout<<"Plug-in: calling gdk_threads_leave()"<<std::endl;
#ifndef _WIN32
  gdk_threads_leave();
#endif
  std::cout<<"Plug-in: gdk_threads_leave() done"<<std::endl;
  return;
}
