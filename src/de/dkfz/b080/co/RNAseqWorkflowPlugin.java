package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 * * TODO Recreate class. Put in dependencies to other workflows, descriptions, capabilities (like ui settings, components) etc.
 */
public class RNAseqWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "2.0.0";
    public static final String CURRENT_VERSION_BUILD_DATE = "Thu Feb 27 14:06:55 CET 2020";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}

