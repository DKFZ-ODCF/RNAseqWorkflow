/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/RNAseqWorkflow).
 */

package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COConfig
import de.dkfz.roddy.core.ExecutionContext
import groovy.transform.CompileStatic

@CompileStatic
class RNAseqConfig extends COConfig {
    RNAseqConfig(ExecutionContext context) {
        super(context)
    }

    boolean usePairedEndProcessing() {
        return !useSingleEndProcessing()
    }

    boolean useSingleEndProcessing() {
        return configValues.getBoolean("useSingleEndProcessing", false)
    }
}
