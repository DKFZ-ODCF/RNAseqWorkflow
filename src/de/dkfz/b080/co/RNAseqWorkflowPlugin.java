package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 * * TODO Recreate class. Put in dependencies to other workflows, descriptions, capabilities (like ui settings, components) etc.
 */
public class RNAseqWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "1.0.11";
    public static final String CURRENT_VERSION_BUILD_DATE = "Mon Jun 12 13:52:13 CEST 2017";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}

