package de.dkfz.b080.co.rnaseqworkflow

import java.lang.reflect.Field

class Helper {

    static void setPrivateField(String name, Object object, Object value) {
        Field f = null
        Class cls = object.class
        while (!f && cls) {
            try {
                f = cls.getDeclaredField(name)
            } catch (Exception ex) {
            }
            cls = cls.superclass
        }
        assert f
        f.setAccessible(true)
        f.set(object, value)
    }

}
