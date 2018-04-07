static void
mask_process(GwyDataField *dfield,
             GwyDataField *existing_mask,
             GwyDataField *maskfield,
             MarkArgs *args)
{
    GwyDataField *output_field;
    gboolean is_field;

    is_field = FALSE;
    output_field = gwy_data_field_new_alike(dfield, FALSE);

    if (args->is_height) {
        gwy_data_field_grains_mark_height(dfield, maskfield, args->height,
                                          args->inverted);
        is_field = TRUE;
    }
    if (args->is_slope) {
        gwy_data_field_grains_mark_slope(dfield, output_field,
                                         args->slope, FALSE);
        if (is_field) {
            if (args->merge_type == GWY_MERGE_UNION)
                gwy_data_field_grains_add(maskfield, output_field);
            else if (args->merge_type == GWY_MERGE_INTERSECTION)
                gwy_data_field_grains_intersect(maskfield, output_field);
        }
        else
            gwy_data_field_copy(output_field, maskfield, FALSE);
        is_field = TRUE;
    }
    if (args->is_lap) {
        gwy_data_field_grains_mark_curvature(dfield, output_field,
                                             args->lap, FALSE);
        if (is_field) {
            if (args->merge_type == GWY_MERGE_UNION)
                gwy_data_field_grains_add(maskfield, output_field);
            else if (args->merge_type == GWY_MERGE_INTERSECTION)
                gwy_data_field_grains_intersect(maskfield, output_field);
        }
        else
            gwy_data_field_copy(output_field, maskfield, FALSE);
    }
    if (existing_mask && args->combine) {
        if (args->combine_type == GWY_MERGE_UNION)
            gwy_data_field_grains_add(maskfield, existing_mask);
        else if (args->combine_type == GWY_MERGE_INTERSECTION)
            gwy_data_field_grains_intersect(maskfield, existing_mask);
    }

    g_object_unref(output_field);
}

/**
 * gwy_data_field_grains_mark_height:
 * @data_field: Data to be used for marking.
 * @grain_field: Data field to store the resulting mask to.
 * @threshval: Relative height threshold, in percents.
 * @below: If %TRUE, data below threshold are marked, otherwise data above
 *         threshold are marked.
 *
 * Marks data that are above/below height threshold.
 **/
void
gwy_data_field_grains_mark_height(GwyDataField *data_field,
                                  GwyDataField *grain_field,
                                  gdouble threshval,
                                  gboolean below)
{
    gdouble min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    gwy_data_field_copy(data_field, grain_field, FALSE);
    gwy_data_field_get_min_max(grain_field, &min, &max);
    if (below)
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 1, 0);
    else
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 0, 1);

    gwy_data_field_invalidate(grain_field);
}

/**
 * gwy_data_field_grains_mark_curvature:
 * @data_field: Data to be used for marking.
 * @grain_field: Data field to store the resulting mask to.
 * @threshval: Relative curvature threshold, in percents.
 * @below: If %TRUE, data below threshold are marked, otherwise data above
 *         threshold are marked.
 *
 * Marks data that are above/below curvature threshold.
 **/
void
gwy_data_field_grains_mark_curvature(GwyDataField *data_field,
                                     GwyDataField *grain_field,
                                     gdouble threshval,
                                     gboolean below)
{
    gdouble min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    gwy_data_field_copy(data_field, grain_field, FALSE);
    gwy_data_field_filter_laplacian(grain_field);

    gwy_data_field_get_min_max(grain_field, &min, &max);
    if (below)
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 1, 0);
    else
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 0, 1);

    gwy_data_field_invalidate(grain_field);
}

/**
 * gwy_data_field_grains_mark_slope:
 * @data_field: Data to be used for marking.
 * @grain_field: Data field to store the resulting mask to.
 * @threshval: Relative slope threshold, in percents.
 * @below: If %TRUE, data below threshold are marked, otherwise data above
 *         threshold are marked.
 *
 * Marks data that are above/below slope threshold.
 **/
void
gwy_data_field_grains_mark_slope(GwyDataField *data_field,
                                 GwyDataField *grain_field,
                                 gdouble threshval,
                                 gboolean below)
{
    GwyDataField *masky;
    gdouble *gdata;
    gint i;
    gdouble xres, yres, min, max;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    g_return_if_fail(GWY_IS_DATA_FIELD(grain_field));

    xres = data_field->xres;
    yres = data_field->yres;

    masky = gwy_data_field_duplicate(data_field);
    gwy_data_field_copy(data_field, grain_field, FALSE);
    gwy_data_field_filter_sobel(grain_field, GWY_ORIENTATION_HORIZONTAL);
    gwy_data_field_filter_sobel(masky, GWY_ORIENTATION_VERTICAL);

    gdata = grain_field->data;
    for (i = 0; i < xres*yres; i++)
        gdata[i] = hypot(gdata[i], masky->data[i]);

    gwy_data_field_get_min_max(grain_field, &min, &max);
    if (below)
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 1, 0);
    else
        gwy_data_field_threshold(grain_field,
                                 min + threshval*(max - min)/100.0, 0, 1);

    g_object_unref(masky);
    gwy_data_field_invalidate(grain_field);
}


/**
 * gwy_data_field_get_min_max:
 * @data_field: A data field.
 * @min: Location to store minimum to.
 * @max: Location to store maximum to.
 *
 * Finds minimum and maximum values of a data field.
 **/
void
gwy_data_field_get_min_max(GwyDataField *data_field,
                           gdouble *min,
                           gdouble *max)
{
    gboolean need_min = FALSE, need_max = FALSE;
    gdouble min1, max1;
    const gdouble *p;
    gint i;

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));

    if (min) {
        if (CTEST(data_field, MIN))
            *min = CVAL(data_field, MIN);
        else
            need_min = TRUE;
    }
    if (max) {
        if (CTEST(data_field, MAX))
            *max = CVAL(data_field, MAX);
        else
            need_max = TRUE;
    }

    if (!need_min && !need_max)
        return;
    else if (!need_min) {
        *max = gwy_data_field_get_max(data_field);
        return;
    }
    else if (!need_max) {
        *min = gwy_data_field_get_min(data_field);
        return;
    }

    min1 = data_field->data[0];
    max1 = data_field->data[0];
    p = data_field->data;
    for (i = data_field->xres * data_field->yres; i; i--, p++) {
        if (G_UNLIKELY(min1 > *p))
            min1 = *p;
        if (G_UNLIKELY(max1 < *p))
            max1 = *p;
    }

    *min = min1;
    *max = max1;
    CVAL(data_field, MIN) = min1;
    CVAL(data_field, MAX) = max1;
    data_field->cached |= CBIT(MIN) | CBIT(MAX);
}


/**
 * gwy_data_field_threshold:
 * @data_field: A data field.
 * @threshval: Threshold value.
 * @bottom: Lower replacement value.
 * @top: Upper replacement value.
 *
 * Tresholds values of a data field.
 *
 * Values smaller than @threshold are set to value @bottom, values higher
 * than @threshold or equal to it are set to value @top
 *
 * Returns: The total number of values above threshold.
 **/
gint
gwy_data_field_threshold(GwyDataField *data_field,
                         gdouble threshval, gdouble bottom, gdouble top)
{
    gint i, n, tot = 0;
    gdouble *p = data_field->data;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field), 0);

    n = data_field->xres * data_field->yres;
    for (i = n; i; i--, p++) {
        if (*p < threshval)
            *p = bottom;
        else {
            *p = top;
            tot++;
        }
    }

    /* We can precompute stats */
    data_field->cached = CBIT(MIN) | CBIT(MAX) | CBIT(SUM) | CBIT(RMS)
                         | CBIT(MED);
    if (tot == n)
        CVAL(data_field, MIN) = CVAL(data_field, MAX) = top;
    else if (tot == 0)
        CVAL(data_field, MIN) = CVAL(data_field, MAX) = bottom;
    else {
        CVAL(data_field, MIN) = MIN(top, bottom);
        CVAL(data_field, MAX) = MAX(top, bottom);
    }
    CVAL(data_field, SUM) = tot*top + (n - tot)*bottom;
    CVAL(data_field, RMS) = (top - bottom)*(top - bottom)
                            * tot/(gdouble)n * (n - tot)/(gdouble)n;
    /* FIXME: may be incorrect for tot == n/2(?) */
    CVAL(data_field, MED) = tot > n/2 ? top : bottom;

    return tot;
}

/**
 * gwy_data_field_filter_laplacian:
 * @data_field: A data field to apply the filter to.
 *
 * Filters a data field with Laplacian filter.
 **/
void
gwy_data_field_filter_laplacian(GwyDataField *data_field)
{
    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_filter_laplacian(data_field, 0, 0,
                                         data_field->xres, data_field->yres);
}

 /**
 * gwy_data_field_area_filter_laplacian:
 * @data_field: A data field to apply the filter to.
 * @col: Upper-left column coordinate.
 * @row: Upper-left row coordinate.
 * @width: Area width (number of columns).
 * @height: Area height (number of rows).
 *
 * Filters a rectangular part of a data field with Laplacian filter.
 **/
void
gwy_data_field_area_filter_laplacian(GwyDataField *data_field,
                                     gint col, gint row,
                                     gint width, gint height)
{
    const gdouble laplace[] = {
        0,  1, 0,
        1, -4, 1,
        0,  1, 0,
    };

    g_return_if_fail(GWY_IS_DATA_FIELD(data_field));
    gwy_data_field_area_convolve_3x3(data_field, laplace,
                                     col, row, width, height);
}

static void
gwy_data_field_area_convolve_3x3(GwyDataField *data_field,
                                 const gdouble *kernel,
                                 gint col, gint row,
                                 gint width, gint height)
{
    gdouble *rm, *rc, *rp;
    gdouble t, v;
    gint xres, i, j;

    xres = data_field->xres;
    rp = data_field->data + row*xres + col;

    /* Special-case width == 1 to avoid complications below.  It's silly but
     * the API guarantees it. */
    if (width == 1) {
        t = rp[0];
        for (i = 0; i < height; i++) {
            rc = rp = data_field->data + (row + i)*xres + col;
            if (i < height-1)
                rp += xres;

            v = (kernel[0] + kernel[1] + kernel[2])*t
                + (kernel[3] + kernel[4] + kernel[5])*rc[0]
                + (kernel[6] + kernel[7] + kernel[8])*rp[0];
            t = rc[0];
            rc[0] = v;
        }
        gwy_data_field_invalidate(data_field);

        return;
    }

    rm = g_new(gdouble, width);
    gwy_assign(rm, rp, width);

    for (i = 0; i < height; i++) {
        rc = rp;
        if (i < height-1)
            rp += xres;
        v = (kernel[0] + kernel[1])*rm[0] + kernel[2]*rm[1]
            + (kernel[3] + kernel[4])*rc[0] + kernel[5]*rc[1]
            + (kernel[6] + kernel[7])*rp[0] + kernel[8]*rp[1];
        t = rc[0];
        rc[0] = v;
        if (i < height-1) {
            for (j = 1; j < width-1; j++) {
                v = kernel[0]*rm[j-1] + kernel[1]*rm[j] + kernel[2]*rm[j+1]
                    + kernel[3]*t + kernel[4]*rc[j] + kernel[5]*rc[j+1]
                    + kernel[6]*rp[j-1] + kernel[7]*rp[j] + kernel[8]*rp[j+1];
                rm[j-1] = t;
                t = rc[j];
                rc[j] = v;
            }
            v = kernel[0]*rm[j-1] + (kernel[1] + kernel[2])*rm[j]
                + kernel[3]*t + (kernel[4] + kernel[5])*rc[j]
                + kernel[6]*rp[j-1] + (kernel[7] + kernel[8])*rp[j];
        }
        else {
            for (j = 1; j < width-1; j++) {
                v = kernel[0]*rm[j-1] + kernel[1]*rm[j] + kernel[2]*rm[j+1]
                    + kernel[3]*t + kernel[4]*rc[j] + kernel[5]*rc[j+1]
                    + kernel[6]*t + kernel[7]*rc[j] + kernel[8]*rc[j+1];
                rm[j-1] = t;
                t = rc[j];
                rc[j] = v;
            }
            v = kernel[0]*rm[j-1] + (kernel[1] + kernel[2])*rm[j]
                + kernel[3]*t + (kernel[4] + kernel[5])*rc[j]
                + kernel[6]*t + (kernel[7] + kernel[8])*rc[j];
        }
        rm[j-1] = t;
        rm[j] = rc[j];
        rc[j] = v;
    }

    g_free(rm);
    gwy_data_field_invalidate(data_field);
}

/**
 * gwy_data_field_grains_add:
 * @grain_field: Field of marked grains (mask).
 * @add_field: Field of marked grains (mask) to be added.
 *
 * Adds @add_field grains to @grain_field.
 *
 * Note: This function is equivalent to
 * |[
 * gwy_data_field_max_of_fields(grain_field, grain_field, add_field);
 * ]|
 **/
void
gwy_data_field_grains_add(GwyDataField *grain_field, GwyDataField *add_field)
{
    gwy_data_field_max_of_fields(grain_field, grain_field, add_field);
}

/**
 * gwy_data_field_max_of_fields:
 * @result: A data field to put the result to.  May be one of @operand1,
 *          @operand2.
 * @operand1: First data field operand.
 * @operand2: Second data field operand.
 *
 * Finds point-wise minima of two data fields.
 **/
void
gwy_data_field_max_of_fields(GwyDataField *result,
                             GwyDataField *operand1,
                             GwyDataField *operand2)
{
    gdouble *p, *q, *r;
    gint xres, yres, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(result));
    g_return_if_fail
        (!gwy_data_field_check_compatibility(result, operand1,
                                            GWY_DATA_COMPATIBILITY_RES));
    g_return_if_fail
        (!gwy_data_field_check_compatibility(result, operand2,
                                            GWY_DATA_COMPATIBILITY_RES));

    xres = result->xres;
    yres = result->yres;
    r = result->data;
    p = operand1->data;
    q = operand2->data;
    for (i = xres*yres; i; i--, p++, q++, r++)
        *r = MAX(*p, *q);

    if (CTEST(operand1, MAX) && CTEST(operand2, MAX)) {
        result->cached = CBIT(MAX);
        CVAL(result, MAX) = MAX(CVAL(operand1, MAX), CVAL(operand2, MAX));
    }
    else
        gwy_data_field_invalidate(result);
}

/**
 * gwy_data_field_grains_intersect:
 * @grain_field: Field of marked grains (mask).
 * @intersect_field: Field of marked grains (mask).
 *
 * Performs intersection betweet two grain fields,
 * result is stored in @grain_field.
 *
 * Note: This function is equivalent to
 * |[
 * gwy_data_field_min_of_fields(grain_field, grain_field, intersect_field);
 * ]|
 **/
void
gwy_data_field_grains_intersect(GwyDataField *grain_field,
                                GwyDataField *intersect_field)
{
    gwy_data_field_min_of_fields(grain_field, grain_field, intersect_field);
}

/**
 * gwy_data_field_min_of_fields:
 * @result: A data field to put the result to.  May be one of @operand1,
 *          @operand2.
 * @operand1: First data field operand.
 * @operand2: Second data field operand.
 *
 * Finds point-wise maxima of two data fields.
 **/
void
gwy_data_field_min_of_fields(GwyDataField *result,
                             GwyDataField *operand1,
                             GwyDataField *operand2)
{
    gdouble *p, *q, *r;
    gint xres, yres, i;

    g_return_if_fail(GWY_IS_DATA_FIELD(result));
    g_return_if_fail
        (!gwy_data_field_check_compatibility(result, operand1,
                                            GWY_DATA_COMPATIBILITY_RES));
    g_return_if_fail
        (!gwy_data_field_check_compatibility(result, operand2,
                                            GWY_DATA_COMPATIBILITY_RES));

    xres = result->xres;
    yres = result->yres;
    r = result->data;
    p = operand1->data;
    q = operand2->data;
    for (i = xres*yres; i; i--, p++, q++, r++)
        *r = MIN(*p, *q);

    if (CTEST(operand1, MIN) && CTEST(operand2, MIN)) {
        result->cached = CBIT(MIN);
        CVAL(result, MIN) = MIN(CVAL(operand1, MIN), CVAL(operand2, MIN));
    }
    else
        gwy_data_field_invalidate(result);
}

/**
 * gwy_data_field_check_compatibility:
 * @data_field1: A data field.
 * @data_field2: Another data field.
 * @check: The compatibility tests to perform.
 *
 * Checks whether two data fields are compatible.
 *
 * Returns: Zero if all tested properties are compatible.  Flags corresponding
 *          to failed tests if data fields are not compatible.
 **/
GwyDataCompatibilityFlags
gwy_data_field_check_compatibility(GwyDataField *data_field1,
                                   GwyDataField *data_field2,
                                   GwyDataCompatibilityFlags check)
{
    GwyDataCompatibilityFlags result = 0;
    gint xres1, xres2, yres1, yres2;
    gdouble xreal1, xreal2, yreal1, yreal2;
    GwySIUnit *unit1, *unit2;

    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field1), check);
    g_return_val_if_fail(GWY_IS_DATA_FIELD(data_field2), check);

    xres1 = data_field1->xres;
    xres2 = data_field2->xres;
    yres1 = data_field1->yres;
    yres2 = data_field2->yres;
    xreal1 = data_field1->xreal;
    xreal2 = data_field2->xreal;
    yreal1 = data_field1->yreal;
    yreal2 = data_field2->yreal;

    /* Resolution */
    if (check & GWY_DATA_COMPATIBILITY_RES) {
        if (xres1 != xres2 || yres1 != yres2)
            result |= GWY_DATA_COMPATIBILITY_RES;
    }

    /* Real size */
    if (check & GWY_DATA_COMPATIBILITY_REAL) {
        /* Keeps the condition in negative form to catch NaNs and odd values
         * as incompatible. */
        if (!(fabs(log(xreal1/xreal2)) <= EPSILON)
            || !(fabs(log(yreal1/yreal2)) <= EPSILON))
            result |= GWY_DATA_COMPATIBILITY_REAL;
    }

    /* Measure */
    if (check & GWY_DATA_COMPATIBILITY_MEASURE) {
        if (!(fabs(log(xreal1/xres1*xres2/xreal2)) <= EPSILON)
            || !(fabs(log(yreal1/yres1*yres2/yreal2)) <= EPSILON))
            result |= GWY_DATA_COMPATIBILITY_MEASURE;
    }

    /* Lateral units */
    if (check & GWY_DATA_COMPATIBILITY_LATERAL) {
        /* This can cause instantiation of data_field units as a side effect */
        unit1 = gwy_data_field_get_si_unit_xy(data_field1);
        unit2 = gwy_data_field_get_si_unit_xy(data_field2);
        if (!gwy_si_unit_equal(unit1, unit2))
            result |= GWY_DATA_COMPATIBILITY_LATERAL;
    }

    /* Value units */
    if (check & GWY_DATA_COMPATIBILITY_VALUE) {
        /* This can cause instantiation of data_field units as a side effect */
        unit1 = gwy_data_field_get_si_unit_z(data_field1);
        unit2 = gwy_data_field_get_si_unit_z(data_field2);
        if (!gwy_si_unit_equal(unit1, unit2))
            result |= GWY_DATA_COMPATIBILITY_VALUE;
    }

    return result;
}