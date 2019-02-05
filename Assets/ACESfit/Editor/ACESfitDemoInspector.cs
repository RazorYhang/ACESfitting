using UnityEngine;
using UnityEditor;
using System.Collections.Generic;

[CustomEditor(typeof(ACESfitDemo))]
public class ACESfitDemoInspector : Editor
{
    public ACESfitDemo demo;
    public override void OnInspectorGUI()
    {
        EditorGUI.BeginChangeCheck();
        DrawDefaultInspector();

        if(EditorGUI.EndChangeCheck())
        {
            if (demo == null)
                demo = this.target as ACESfitDemo;

            demo.debugFit = true;
        }
    }
}
